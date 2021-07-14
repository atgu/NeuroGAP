# This script was written by Zan Koenig for the purpose of conducting Quality Control on the pilot dataset for the
# NeuroGAP-Psychosis project.
# If planning to run Autosomal and X QC, this script is intended to be run second.
# Please see the readme for information on input format and requirements


import hail as hl
from hail.plot import show
hl.plot.output_notebook()
hl.init()


'''
Annotating Functions
'''


# Conducts sex imputation statistics for a site. Returns an mt with an annotated
# imputed sex column & a column which flags those who failed sex filter as True
def check_sex(mt):
    '''
    Conducts sex imputation statistics for a site. Returns an mt with an annotated
    imputed sex column & a column which flags those who failed sex filter as True
    :param mt: hail matrix table which contains a reported sex column of 'F', 'M', or 'U' named "reported_sex"
    :return: hail matrix table with a new column named sex_filter containing the sex discrepancy filter flag
    '''
    new_mt = hl.impute_sex(mt.GT)
    mt = mt.annotate_cols(imputedSex=new_mt[mt.s])
    return mt.annotate_cols(sex_filter=mt.imputedSex.is_female != (mt.reported_sex == 'F'))


# Annotates the mt with the results of the sample call rate filter
# Creates a column which flags those who failed the filter as True
def check_sample_cr_filter(mt, sample_geno_filter):
    '''
    Annotates the mt with the results of the sample call rate filter
    Creates a column which flags those who failed the filter as True
    :param mt: hail matrix table which has run the hl.sample_qc() function with the name param as 'sample_qc'
    :param sample_geno_filter: genotype filter percentage cutoff as a float eg. .98
    :return: hail matrix table with a new column named sample_cr_filter containing the sample call rate filter flag
    '''
    return mt.annotate_cols(sample_cr_filter=mt.sample_qc.call_rate < sample_geno_filter)


# Annotates a matrix table with the following meta data: collaborator participant id, site id, and
# autocall call rate. Assumes that meta and mt don't have the same key, keys meta by chip_well_barcode
def annotate_meta(mt, meta, key_by):
    '''
    Annotates a matrix table with the following meta data: collaborator participant id, site id, and
    autocall call rate. Assumes that meta and mt don't have the same key, keys meta by chip_well_barcode
    :param mt: hail matrix table
    :param meta_table: tsv file containing sample IDs, reported sex information, and IS specific to
    :param key_by:
    :return:
    '''
    # Setting the key of metaDat as chip_well_barcode
    meta = meta.key_by(key_by)

    # Adding the reported sex column from the sample data table to the matrix table
    mt = mt.annotate_cols(reported_sex=meta[mt.s].reported_gender)

    # Annotating the matrix table with the projID from meta data
    mt = mt.annotate_cols(collab_PID=meta[mt.s].collaborator_participant_id)

    # Creating a new column in mt for siteID which is the 3 letters from the collab participant id
    mt = mt.annotate_cols(siteID=mt.col.collab_PID[0:3])

    # Switching lowercase letters to capital in collabPID and siteID
    mt = mt.annotate_cols(collab_PID=mt.collab_PID.upper())
    mt = mt.annotate_cols(siteID=mt.siteID.upper())

    # Annotating mt with autocall call rate 
    mt = mt.annotate_cols(autocall_call_rate=meta[mt.s].autocall_call_rate)

    # Filtering out individuals with an autocall call rate below .95
    mt = mt.filter_cols(hl.float(mt.autocall_call_rate) < .95, keep=False)

    # Filtering NA12878 individual from the matrix table
    samples_to_remove = {'NA12878'}
    set_to_remove = hl.literal(samples_to_remove)
    return mt.filter_cols(~set_to_remove.contains(mt['collab_PID']))


'''
Printing Functions
'''


def print_count(mt):
    '''
    Prints out total sample/variant count for an mt
    :param mt: hail matrix table
    :return: print statement with number of samples and variants
    '''
    n = mt.count()
    print('Number of Samples: {}\nNumber of Variants: {}'.format(n[1], n[0]))


def print_filter_counts(filter_name, mt, count, filter_type):
    '''
    Prints how many remaining samples/variants there are after filtering.
    :param filter_name: name of the filter for which you'd like the pre/post counts
    :param mt: filtered hail matrix table
    :param count: mt.count() pre-filter
    :param filter_type: "samples" for sample filters, "variants" for variant filters
    :return: print statement with the number of samples or variants before & after the specified filter
    '''
    if filter_type == 'variants':
        print('After %s filter, %d/%d %s remain.\n' % (filter_name, mt.count_rows(), count, filter_type))
    elif filter_type == 'samples':
        print('After %s filter, %d/%d %s remain.\n' % (filter_name, mt.count_cols(), count, filter_type))



def site_var_count(mt):
    '''
    Prints out the variant counts per site (hardcoded for NGAP dataset site locations)
    :param mt: NGAP hail matrix table with info from 5 pilot data sites
    :return: print statement with variant counts for each site
    '''
    row_count = {'MOP': 0, 'MAP': 0, 'KWP': 0, 'AAP': 0, 'CTP': 0}
    for site, count in row_count.items():
        mt_tmp = mt.filter_cols(mt.siteID == site)
        row_count[site] = mt_tmp.count_rows()
    print(('Site Variant Counts:\n\tKenya Moi - {}\n\tUganda - {}\n\tKenya Kemri - {}\n\tEthiopia - {}\n\tSouth '
           'Africa - {}').format(row_count['MOP'], row_count['MAP'], row_count['KWP'], row_count['AAP'],
                                 row_count['CTP']))


# Prints out the sample counts for each site
def site_sample_count(mt):
    '''
    Prints out the sample counts per site (hardcoded for NGAP dataset site locations)
    :param mt: NGAP hail matrix table with info from 5 pilot data sites
    :return: print statement with sample counts for each site
    '''
    count = mt.aggregate_cols(hl.agg.counter(mt.siteID))
    print('Site Sample Counts:')
    for key in count.keys():
        if key == "MOP":
            print('\tKenya Moi - {}'.format(count["MOP"]))
        elif key == "MAP":
            print('\tUganda - {}'.format(count["MAP"]))
        elif key == 'KWP':
            print('\tKenya Kemri - {}'.format(count["KWP"]))
        elif key == "AAP":
            print('\tEthiopia - {}'.format(count["AAP"]))
        elif key == "CTP":
            print('\tSouth Africa - {}'.format(count["CTP"]))


def print_site_indel_snp(mt):
    '''
    Prints out the average indel and SNP counts for each site
    :param mt: NGAP hail matrix table with info from 5 pilot data sites and hl.sample_qc() run
    :return: print statement with average indel and SNP counts for each site
    '''
    site_ids = ['MOP', 'MAP', 'KWP', 'AAP', 'CTP']
    for i in site_ids:
        site_mt = mt.filter_cols(mt.siteID == i)
        del_avg = site_mt.aggregate_cols(hl.agg.sum(site_mt.sample_qc.n_deletion)) / site_mt.count_cols()
        ins_avg = site_mt.aggregate_cols(hl.agg.sum(site_mt.sample_qc.n_insertion)) / site_mt.count_cols()
        snp_avg = site_mt.aggregate_cols(hl.agg.sum(site_mt.sample_qc.n_snp)) / site_mt.count_cols()
        if i == "MOP":
            print(
                ('Kenya Moi Average Indel and SNP Count:\n\tAverage Indel Count - {:.0f}\n\tAverage SNP count - {:.0f}'
                 '\n').format(
                    (del_avg + ins_avg), snp_avg))
        elif i == "MAP":
            print(
                ('Uganda Average Indel and SNP Count:\n\tAverage Indel Count - {:.0f}\n\tAverage SNP count - {:.0f}'
                 '\n').format(
                    (del_avg + ins_avg), snp_avg))
        elif i == 'KWP':
            print(
                ('Kenya Kemri Average Indel and SNP Count:\n\tAverage Indel Count - {:.0f}\n\tAverage SNP count - '
                 '{:.0f}''\n').format(
                    (del_avg + ins_avg), snp_avg))
        elif i == "AAP":
            print(
                ('Ethiopia Average Indel and SNP Count:\n\tAverage Indel Count - {:.0f}\n\tAverage SNP count - {:.0f}'
                 '\n').format(
                    (del_avg + ins_avg), snp_avg))
        elif i == "CTP":
            print(
                ('South Africa Average Indel and SNP Count:\n\tAverage Indel Count - {:.0f}\n\tAverage SNP count - '
                 '{:.0f}''\n').format((del_avg + ins_avg), snp_avg))


def print_indel_snp(mt):
    '''
    Prints out the average indel and SNP counts for all sites combined
    :param mt: hail matrix table with hl.sample_qc() run
    :return: print statement with average indel and SNP counts for the dataset
    '''
    del_avg = mt.aggregate_cols(hl.agg.sum(mt.sample_qc.n_deletion)) / mt.count_cols()
    ins_avg = mt.aggregate_cols(hl.agg.sum(mt.sample_qc.n_insertion)) / mt.count_cols()
    snp_avg = mt.aggregate_cols(hl.agg.sum(mt.sample_qc.n_snp)) / mt.count_cols()
    print(
        ('Average Indel and SNP Counts for all Sites:\n\tAverage Indel Count - {:.0f}\n\tAverage SNP count - {:.0f}'
         '\n').format((del_avg + ins_avg), snp_avg))


'''
Formatting Functions
'''


def align_mt2_to_mt1(mt1, mt2):
    '''
    Takes on hail matrix table and aligns to another based on col_index
    :param mt1: hail matrix table with desired column index
    :param mt2: hail matrix table to be aligned given mt1 col index
    :return: hail matrix table with data from mt2 rearranged to the col order of mt1
    '''
    mt2 = mt2.add_col_index()
    new_col_order = mt2.index_cols(mt1.col_key).col_idx.collect()
    return mt2.choose_cols(new_col_order)


# Path for metadata file containing reported sex, collaborator participant id
meta_data = ('gs://neurogap/Pilot_Data_HailQC/sample_AddisEthiopia_KEMRIKenya_MakerereUganda_MoiKenya_UCTSouthAfrica_.'
             'tsv')

# List of all the paths to the vcfs for each geographic location (site)
mt_paths = [
    ('gs://fc-087888c3-eb94-4828-9926-d8d713eef821/NeuroGap_Koenen_Atwoli_Kenya-MoiPilot_GSA-MD_v1/NeuroGap_Koenen_'
     'Atwoli_Kenya-MoiPilot_GSA-MD_v1.vcf.gz'),
    ('gs://fc-cd517a08-fbb2-451f-a104-d21b32c3663c/NeuroGap_Koenen_Akena_UgandaPilot_GSA-MD_v1/NeuroGap_Koenen_Akena_'
     'UgandaPilot_GSA-MD_v1.vcf.gz'),
    ('gs://fc-abcbfa13-e1ff-4edb-bbb9-8ad757b3f2e0/NeuroGap_Koenen_Newton__Kenya-KEMRI_GSA-MD/NeuroGap_Koenen_Newton__'
     'Kenya-KEMRI_GSA-MD.vcf.gz'),
    ('gs://fc-8b77bd6a-19a6-42aa-9670-867fbd768b74/NeuroGap_Koenen_Teferra_Ethiopia_GSA-MD/NeuroGap_Koenen_Teferra_'
     'Ethiopia_GSA-MD.vcf.gz'),
    ('gs://fc-b2479a12-eaf2-4f96-891f-920a0191c252/NeuroGap_Koenen_Stein_SouthAfrica_GSA-MD/NeuroGap_Koenen_Stein_'
     'SouthAfrica_GSA-MD.vcf.gz')]


# Reading in and creating a list of all of the site matrix tables 
mt_list = [hl.import_vcf(mt_path, force_bgz=True) for mt_path in mt_paths]

# Importing the metadata file as a hail table
meta = hl.import_table(meta_data)

# Reading in the postQC matrix table and subsetting to autosomal data
mt_auto = hl.read_matrix_table('gs://neurogap/Pilot_Data_HailQC/NeuroGAP_pilotAuto_clean.mt')
mt_preQC = hl.read_matrix_table('gs://neurogap/Pilot_Data_HailQC/pilotMT_preQC.mt')


# Annotating male matrix tables with metadata from the meta table (see annotate_meta() for details)
mt_list = [annotate_meta(mt, meta, 'chip_well_barcode') for mt in mt_list]

# Filtering to X raw chromosome data
mt_x_list = [mt.filter_rows(mt.locus.contig == 'X') for mt in mt_list]

# From the X data, splitting by sex
female_list = [mt.filter_cols(mt.reported_sex == 'F') for mt in mt_x_list]
male_list = [mt.filter_cols(mt.reported_sex == 'M') for mt in mt_x_list]

# Annotating the female matrix tables with variant QC data
female_list = [hl.variant_qc(mt, name='variant_qc') for mt in female_list]

# Annotating the male matrix tables with variant QC data
male_list = [hl.variant_qc(mt, name='variant_qc') for mt in male_list]

# Annotating the female matrix tables with -log10 hwe pval data
female_list = [mt.annotate_rows(log10_hwe_pval=-hl.log10(mt.variant_qc.p_value_hwe)) for mt in female_list]

# Annotating the male matrix tables with -log10 hwe pval data
male_list = [mt.annotate_rows(log10_hwe_pval=-hl.log10(mt.variant_qc.p_value_hwe)) for mt in male_list]


'''
Plotting female hwe pval distributions to understand X chrom QC

 MOP = Kenya Moi
 AAP = Ethiopia
 KWP = Kenya Kemri
 CTP = South Africa
 MAP = Uganda
'''

# Plotting Kenya Moi
p = hl.plot.histogram(female_list[0].log10_hwe_pval, legend="HWE -log10 p-val", range=(0, 20),
                      title="Kenya, Moi -log10 HWE p-val Before Prelim QC")
show(p)

# Plotting Ethiopia
p = hl.plot.histogram(female_list[1].log10_hwe_pval, legend="HWE -log10 p-val", range=(0, 20),
                      title="Uganda -log10 HWE p-val Before Prelim QC")
show(p)

# Plotting KEMRI
p = hl.plot.histogram(female_list[2].log10_hwe_pval, legend="HWE -log10 p-val", range=(0, 20),
                      title="KEMRI -log10 HWE p-val Before Prelim QC")
show(p)

# Plotting South Africa
p = hl.plot.histogram(female_list[3].log10_hwe_pval, legend="HWE -log10 p-val", range=(0, 20),
                      title="Ethiopia -log10 HWE p-val Before Prelim QC")
show(p)

# Plotting Uganda
p = hl.plot.histogram(female_list[4].log10_hwe_pval, legend="HWE -log10 p-val", range=(0, 20),
                      title="South Africa -log10 HWE p-val Before Prelim QC")
show(p)


# Joining the male matrix tables using union_cols(), will need to then to annotate with row data 
for i in range(len(male_list) - 1):
    if i == 0:
        mt_m = male_list[i].union_cols(male_list[i + 1])
    else:
        mt0 = mt_m.union_cols(male_list[i + 1])
        mt_m = mt0


# Joining the female matrix tables using union_cols(), will need to then to annotate with row data 
for i in range(len(mt_list) - 1):
    if i == 0:
        mt_f = female_list[i].union_cols(female_list[i + 1])
    else:
        mt0 = mt_f.union_cols(female_list[i + 1])
        mt_f = mt0


# Filtering the raw X chrom data to just individuals who passed QC
auto_cols = mt_auto.cols()
mt_m = mt_m.semi_join_cols(auto_cols)


# Filtering the raw X chrom data to just individuals who passed QC
auto_cols = mt_auto.cols()
mt_f = mt_f.semi_join_cols(auto_cols)


mt_f.write('gs://neurogap/Pilot_Data_HailQC/NeuroGAP_female_x.mt')


mt_m.write('gs://neurogap/Pilot_Data_HailQC/NeuroGAP_male_x.mt')

mt_f = hl.read_matrix_table('gs://neurogap/Pilot_Data_HailQC/NeuroGAP_female_x.mt')
mt_m = hl.read_matrix_table('gs://neurogap/Pilot_Data_HailQC/NeuroGAP_male_x.mt')


# Creating a list of site IDs to annotate globals 
site_ids = hl.array(['MOP', 'MAP', 'KWP', 'AAP ', 'CTP'])

mt_f = mt_f.annotate_globals(location=site_ids)
mt_m = mt_m.annotate_globals(location=site_ids)


# Filtering the females to par and non par
intervals = [hl.parse_locus_interval(x) for x in ['X:60001-2699520', 'X:154931044-155260560']]
femaleX_noPar = hl.filter_intervals(mt_f, intervals, keep=False)
femaleX_par = hl.filter_intervals(mt_f, intervals, keep=True)


# Filtering the males to par and non par
intervals = [hl.parse_locus_interval(x) for x in ['X:60001-2699520', 'X:154931044-155260560']]
maleX_noPar = hl.filter_intervals(mt_m, intervals, keep=False)
maleX_par = hl.filter_intervals(mt_m, intervals, keep=True)


'''
Conducting variant QC for the PAR region of the X chrom
'''

# Merging the female PAR region with the male PAR region
mt_x_par = femaleX_par.union_cols(maleX_par)


# Annotating the matrix tables with variant QC data
mt_x_list = [hl.variant_qc(mt, name='variant_qc') for mt in mt_x_list]


intervals = [hl.parse_locus_interval(x) for x in ['X:60001-2699520', 'X:154931044-155260560']]
mt_x_list = [hl.filter_intervals(mt, intervals, keep=True) for mt in mt_x_list]

# Creating lists for variant call rate, hwe pval and maf in joined mt
mt_x_par = mt_x_par.annotate_rows(var_call_rate=hl.empty_array('float64'))
mt_x_par = mt_x_par.annotate_rows(hwe_pval=hl.empty_array('float64'))
mt_x_par = mt_x_par.annotate_rows(maf=hl.empty_array('float64'))


# Annotating var_call_rate list with variant call rates
for mt_next in mt_x_list:
    mt_x_par = mt_x_par.annotate_rows(
        var_call_rate=mt_x_par.var_call_rate.append(mt_next.index_rows(mt_x_par.row_key).variant_qc.call_rate))

# Annotating hwe_pval list with hwe pvals
for mt_next in mt_x_list:
    mt_x_par = mt_x_par.annotate_rows(
        hwe_pval=mt_x_par.hwe_pval.append((mt_next.index_rows(mt_x_par.row_key).variant_qc.p_value_hwe)))

# Annotating maf list with mafs
for mt_next in mt_x_list:
    mt_x_par = mt_x_par.annotate_rows(
        maf=mt_x_par.maf.append(hl.min(mt_next.index_rows(mt_x_par.row_key).variant_qc.AF)))


# Annotating joined mt_x_par with lists for var call rate, hwe pval and maf flags
mt_x_par = mt_x_par.annotate_rows(var_cr_flag=hl.empty_array('bool'))
mt_x_par = mt_x_par.annotate_rows(hwe_pval_flag=hl.empty_array('bool'))
mt_x_par = mt_x_par.annotate_rows(maf_flag=hl.empty_array('bool'))


# Filling var call rate flag list
for mt_next in mt_x_list:
    mt_x_par = mt_x_par.annotate_rows(
        var_cr_flag=mt_x_par.var_cr_flag.append(mt_next.index_rows(mt_x_par.row_key).variant_qc.call_rate < .95))

# Filling hwe pval flag list
for mt_next in mt_x_list:
    mt_x_par = mt_x_par.annotate_rows(hwe_pval_flag=mt_x_par.hwe_pval_flag.append(
        mt_next.index_rows(mt_x_par.row_key).variant_qc.p_value_hwe <= 1e-03))

# Filling maf flag list
for mt_next in mt_x_list:
    mt_x_par = mt_x_par.annotate_rows(
        maf_flag=mt_x_par.maf_flag.append(hl.min(mt_next.index_rows(mt_x_par.row_key).variant_qc.AF) <= 0.005))


'''
Running QC on the PAR region of the X Chromosome
'''


# SNP call rate 1st pass filtering 
par_qc = mt_x_par.filter_rows((mt_x_par.var_cr_flag.contains(True) == True), keep=False)

# Printing out counts post filter
print_filter_counts('SNP call rate', par_qc, mt_x_par.count_rows(), 'variants')


# Variable for count pre filter
count = par_qc.count_rows()

# MAF filtering
par_qc = par_qc.filter_rows((par_qc.maf_flag.contains(True) == True), keep=False)

# # Printing out counts post filter
print_filter_counts('MAF', par_qc, count, 'variants')


# Variable for count pre filter
count = par_qc.count_rows()

# HWE filtering
par_qc = par_qc.filter_rows((par_qc.hwe_pval_flag.contains(True) == True), keep=False)

# Printing out counts post filter
print_filter_counts('HWE', par_qc, count, 'variants')

par_qc.write('gs://neurogap/Pilot_Data_HailQC/NeuroGAP_PAR_clean.mt')


'''
Adding filter info for female non-par matrix table
'''


# Creating lists for variant call rate, hwe pval and maf in joined mt 
femaleX_noPar = femaleX_noPar.annotate_rows(var_call_rate=hl.empty_array('float64'))
femaleX_noPar = femaleX_noPar.annotate_rows(hwe_pval=hl.empty_array('float64'))
femaleX_noPar = femaleX_noPar.annotate_rows(maf=hl.empty_array('float64'))


# Annotating var_call_rate list with variant call rates
for mt_next in female_list:
    femaleX_noPar = femaleX_noPar.annotate_rows(var_call_rate=femaleX_noPar.var_call_rate.append(
        mt_next.index_rows(femaleX_noPar.row_key).variant_qc.call_rate))


# Annotating hwe_pval list with hwe pvals
for mt_next in female_list:
    femaleX_noPar = femaleX_noPar.annotate_rows(
        hwe_pval=femaleX_noPar.hwe_pval.append((mt_next.index_rows(femaleX_noPar.row_key).variant_qc.p_value_hwe)))


# Annotating maf list with mafs
for mt_next in female_list:
    femaleX_noPar = femaleX_noPar.annotate_rows(
        maf=femaleX_noPar.maf.append(hl.min(mt_next.index_rows(femaleX_noPar.row_key).variant_qc.AF)))


# Annotating joined mt with lists for var call rate, hwe pval and maf flags
femaleX_noPar = femaleX_noPar.annotate_rows(var_cr_flag=hl.empty_array('bool'))
femaleX_noPar = femaleX_noPar.annotate_rows(hwe_pval_flag=hl.empty_array('bool'))
femaleX_noPar = femaleX_noPar.annotate_rows(maf_flag=hl.empty_array('bool'))


# Filling var call rate flag list for X chrom QC
for mt_next in female_list:
    femaleX_noPar = femaleX_noPar.annotate_rows(var_cr_flag=femaleX_noPar.var_cr_flag.append(
        (mt_next.index_rows(femaleX_noPar.row_key).variant_qc.call_rate) < .98))


# Filling maf flag list for X chrom QC
for mt_next in female_list:
    femaleX_noPar = femaleX_noPar.annotate_rows(
        maf_flag=femaleX_noPar.maf_flag.append(hl.min(mt_next.index_rows(femaleX_noPar.row_key).variant_qc.AF) <= 0.01))


# Filling hwe pval flag list for X chrom QC
for mt_next in female_list:
    femaleX_noPar = femaleX_noPar.annotate_rows(hwe_pval_flag=femaleX_noPar.hwe_pval_flag.append(
        (mt_next.index_rows(femaleX_noPar.row_key).variant_qc.p_value_hwe) <= 1e-06))


'''
Adding filter info for male non-par matrix table
'''


# Creating lists for variant call rate, hwe pval and maf in joined mt 
mt_m = mt_m.annotate_rows(var_call_rate=hl.empty_array('float64'))
mt_m = mt_m.annotate_rows(hwe_pval=hl.empty_array('float64'))
mt_m = mt_m.annotate_rows(maf=hl.empty_array('float64'))

# Annotating var_call_rate list with variant call rates
for mt_next in male_list:
    mt_m = mt_m.annotate_rows(
        var_call_rate=mt_m.var_call_rate.append(mt_next.index_rows(mt_m.row_key).variant_qc.call_rate))


# Annotating hwe_pval list with hwe pvals
for mt_next in male_list:
    mt_m = mt_m.annotate_rows(hwe_pval=mt_m.hwe_pval.append((mt_next.index_rows(mt_m.row_key).variant_qc.p_value_hwe)))


# Annotating maf list with mafs
for mt_next in male_list:
    mt_m = mt_m.annotate_rows(maf=mt_m.maf.append(hl.min(mt_next.index_rows(mt_m.row_key).variant_qc.AF)))


# Annotating joined mt with lists for var call rate, hwe pval and maf flags
mt_m = mt_m.annotate_rows(var_cr_flag=hl.empty_array('bool'))
mt_m = mt_m.annotate_rows(hwe_pval_flag=hl.empty_array('bool'))
mt_m = mt_m.annotate_rows(maf_flag=hl.empty_array('bool'))


# Filling var call rate flag list for X chrom QC
for mt_next in male_list:
    mt_m = mt_m.annotate_rows(
        var_cr_flag=mt_m.var_cr_flag.append((mt_next.index_rows(mt_m.row_key).variant_qc.call_rate) < .98))


# Filling maf flag list for X chrom QC
for mt_next in male_list:
    mt_m = mt_m.annotate_rows(
        maf_flag=mt_m.maf_flag.append(hl.min(mt_next.index_rows(mt_m.row_key).variant_qc.AF) <= 0.01))


# Filling hwe pval flag list for X chrom QC
for mt_next in male_list:
    mt_m = mt_m.annotate_rows(
        hwe_pval_flag=mt_m.hwe_pval_flag.append((mt_next.index_rows(mt_m.row_key).variant_qc.p_value_hwe) <= 1e-06))


'''
Female X Chromosome no Par Filtering
'''

# X Chrom SNP call rate filtering 
fx_noPar_qc = femaleX_noPar.filter_rows((femaleX_noPar.var_cr_flag.contains(True) == True), keep=False)

# Printing out counts post filter
print_filter_counts('SNP call rate', fx_noPar_qc, mt_f.count_rows(), 'variants')


# Variable for count pre filter
count = fx_noPar_qc.count_rows()

# X Chrom MAF filtering
fx_noPar_qc = fx_noPar_qc.filter_rows((fx_noPar_qc.maf_flag.contains(True) == True), keep=False)

# Printing out counts post filter
print_filter_counts('MAF', fx_noPar_qc, count, 'variants')


# Variable for count pre filter
count = fx_noPar_qc.count_rows()

# HWE filtering
fx_noPar_qc = fx_noPar_qc.filter_rows((fx_noPar_qc.hwe_pval_flag.contains(True) == True), keep=False)

# Printing out counts post filter
print_filter_counts('HWE', fx_noPar_qc, count, 'variants')


'''
Male X nonPar Chromosome Filtering
'''


# X Chrom SNP call rate filtering 
mt_qc_m = mt_m.filter_rows((mt_m.var_cr_flag.contains(True) == True), keep=False)

# Printing out counts post filter
print_filter_counts('SNP call rate', mt_qc_m, mt_m.count_rows(), 'variants')

# Variable for count pre filter
count = mt_qc_m.count_rows()

# X Chrom MAF filtering
mt_qc_m = mt_qc_m.filter_rows((mt_qc_m.maf_flag.contains(True) == True), keep=False)

# Printing out counts post filter
print_filter_counts('MAF', mt_qc_m, count, 'variants')


# Variable for count pre filter
count = mt_qc_m.count_rows()

# HWE filtering
mt_qc_m = mt_qc_m.filter_rows((mt_qc_m.hwe_pval_flag.contains(True) == True), keep=False)

# Printing out counts post filter
print_filter_counts('HWE', mt_qc_m, count, 'variants')

# Writing out female and male X chromosome mts post QC
mt_qc_m.write('gs://neurogap/Pilot_Data_HailQC/NeuroGAP_MaleX_clean.mt')

fx_noPar_qc.write('gs://neurogap/Pilot_Data_HailQC/NeuroGAP_FemaleXnonPAR_clean.mt')

# Reading in all the necessary mts to combine into final dataset
mt_auto = hl.read_matrix_table('gs://neurogap/Pilot_Data_HailQC/NeuroGAP_pilotAuto_clean.mt')

fnonPAR_qc = hl.read_matrix_table('gs://neurogap/Pilot_Data_HailQC/NeuroGAP_FemaleXnonPAR_clean.mt')

mt_preQC = hl.read_matrix_table('gs://neurogap/Pilot_Data_HailQC/pilotMT_preQC.mt')

mt_par = hl.read_matrix_table('gs://neurogap/Pilot_Data_HailQC/NeuroGAP_PAR_clean.mt')

mt_x_m = hl.read_matrix_table('gs://neurogap/Pilot_Data_HailQC/NeuroGAP_male_x.mt')


# Filtering the males to par and non par
intervals = [hl.parse_locus_interval(x) for x in ['X:60001-2699520', 'X:154931044-155260560']]
maleX_noPar = hl.filter_intervals(mt_x_m, intervals, keep=False)


print_count(maleX_noPar)


# Adding male nonpar to the female nonpar post qc 
nonpar_f_rows = fnonPAR_qc.rows()
nonpar_m = maleX_noPar.semi_join_rows(nonpar_f_rows)


print_count(nonpar)


nonpar = nonpar_m.union_cols(fnonPAR_qc)


# Sanity checking the var/sample counts for each dataset
print_count(nonpar)
print_count(mt_preQC)
print_count(mt_auto)
print_count(fnonPAR_qc)
print_count(mt_par)


# Changing the keys for all the mts to match
fnonPAR_qc = fnonPAR_qc.key_rows_by(fnonPAR_qc.rsid)
mt_par = mt_par.key_rows_by(mt_par.rsid)
mt_auto = mt_auto.key_rows_by(mt_auto.rsid)
nonpar = nonpar.key_rows_by(nonpar.rsid)


# Dropping fields added during qc from the datasets pre-merge
nonpar = nonpar.drop('log10_hwe_pval')
mt_par = mt_par.drop('log10_hwe_pval')
mt_auto = mt_auto.drop('sample_qc', 'imputedSex', 'sex_filter', 'sample_cr_filter',
                       'related_filter')

# Aligning the column order of mt_par to mt_auto
mt_par = align_mt2_to_mt1(mt_auto, mt_par)

# Joining mt_auto and mt_par
mt_clean = mt_auto.union_rows(mt_par)

# Sanity check counts
print_count(mt_clean)

# Dropping unnecessary field from mt
fnonPAR_qc = fnonPAR_qc.drop('log10_hwe_pval')

# Aligning mt_clean col order to nonpar
mt_clean = align_mt2_to_mt1(nonpar, mt_clean)

# Dropping fields added during qc from mt
mt_clean = mt_clean.drop('var_call_rate', 'hwe_pval', 'maf',
                         'var_cr_flag', 'hwe_pval_flag', 'maf_flag')

# Joining mt_clean and nonpar
mt_clean = nonpar.union_rows(mt_clean)

# Sanity check counts
print_count(mt_clean)

# Changing sex annotation so it will properly be output when converted to plink format
mt_clean = mt_clean.annotate_cols(sex=hl.if_else((mt_clean.reported_sex == 'F'), True, False))

# Sanity check sex count after changing format
mt_clean.aggregate_cols(hl.agg.counter(mt_clean.sex == True))

# Reading in csv with FID info
fid = hl.import_table('gs://neurogap/Pilot_Data_HailQC/fid_info.csv', delimiter=',')

# Changing keys to match mt_clean
fid = fid.key_by(fid.siteID)

# Adding on proper FIDs to dataset (for plink output)
mt_clean = mt_clean.annotate_cols(fid=fid[mt_clean.siteID].FID)

# Changing row key to what we want for plink output
mt_clean = mt_clean.key_rows_by(locus=mt_clean['locus'], alleles=mt_clean['alleles'])

# Output to plink, specifying desired fields
hl.export_plink(mt_clean, 'gs://neurogap-pilot-clean/NeuroGAP_pilot_clean', ind_id=mt_clean.collab_PID,
                fam_id=mt_clean.fid, is_female=mt_clean.sex, varid=mt_clean.rsid)
