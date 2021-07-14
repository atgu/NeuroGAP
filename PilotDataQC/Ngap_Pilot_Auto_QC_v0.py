# This script was written by Zan Koenig for the purpose of conducting Quality Control on the pilot dataset for the
# NeuroGAP-Psychosis project.
# If planning to run Autosomal and X QC, this script is intended to be run first.
# Please see the readme for information on input format and requirements


import hail as hl
hl.plot.output_notebook()
hl.init()


'''
Annotating Functions
'''


def check_sex(mt):
    '''
    Conducts sex imputation statistics for a site. Returns an mt with an annotated
    imputed sex column & a column which flags those who failed sex filter as True
    :param mt: hail matrix table which contains a reported sex column of 'F', 'M', or 'U' named "reported_sex"
    :return: hail matrix table with a new column named sex_filter containing the sex discrepancy filter flag
    '''
    new_mt = hl.impute_sex(mt.GT)
    mt = mt.annotate_cols(imputedSex=new_mt[mt.s])
    return (
        mt.annotate_cols(sex_filter=((mt.imputedSex.is_female != (mt.reported_sex == 'F')) | (mt.reported_sex == 'U'))))


def check_sample_cr_filter(mt, sample_geno_filter):
    '''
    Annotates the mt with the results of the sample call rate filter
    Creates a column which flags those who failed the filter as True
    :param mt: hail matrix table which has run the hl.sample_qc() function with the name param as 'sample_qc'
    :param sample_geno_filter: genotype filter percentage cutoff as a float eg. .98
    :return: hail matrix table with a new column named sample_cr_filter containing the sample call rate filter flag
    '''
    return mt.annotate_cols(sample_cr_filter=mt.sample_qc.call_rate < sample_geno_filter)


def annotate_meta(mt, meta_table, key_by):
    '''
    Annotates a matrix table with the following meta data: collaborator participant id, site id, and
    autocall call rate. Assumes that meta and mt don't have the same key, keys meta by chip_well_barcode
    :param mt: hail matrix table
    :param meta_table: tsv file containing sample IDs, reported sex information, and IS specific to
    :param key_by:
    :return:
    '''
    # Setting the key of metaDat as chip_well_barcode
    meta_table = meta_table.key_by(key_by)

    # Adding the reported sex column from the sample data table to the matrix table
    mt = mt.annotate_cols(reported_sex=meta_table[mt.s].reported_gender)

    # Annotating the matrix table with the projID from meta data
    mt = mt.annotate_cols(collab_PID=meta_table[mt.s].collaborator_participant_id)

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


def print_filter_counts(filter_name, mt, count_num, filter_type):
    '''
    Prints how many remaining samples/variants there are after filtering.
    :param filter_name: name of the filter for which you'd like the pre/post counts
    :param mt: filtered hail matrix table
    :param count: mt.count() pre-filter
    :param filter_type: "samples" for sample filters, "variants" for variant filters
    :return: print statement with the number of samples or variants before & after the specified filter
    '''
    if filter_type == 'variants':
        print('After %s filter, %d/%d %s remain.\n' % (filter_name, mt.count_rows(), count_num, filter_type))
    elif filter_type == 'samples':
        print('After %s filter, %d/%d %s remain.\n' % (filter_name, mt.count_cols(), count_num, filter_type))


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
    print(
        ('Site Variant Counts:\n\tKenya Moi - {}\n\tUganda - {}\n\tKenya Kemri - {}\n\tEthiopia - {}\n\tSouth Africa '
         '- {}').format(
            row_count['MOP'], row_count['MAP'], row_count['KWP'], row_count['AAP'], row_count['CTP']))


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
    siteIDs = ['MOP', 'MAP', 'KWP', 'AAP', 'CTP']
    for i in siteIDs:
        site_mt = mt.filter_cols(mt.siteID == i)
        del_avg = site_mt.aggregate_cols(hl.agg.sum(site_mt.sample_qc.n_deletion)) / site_mt.count_cols()
        ins_avg = site_mt.aggregate_cols(hl.agg.sum(site_mt.sample_qc.n_insertion)) / site_mt.count_cols()
        snp_avg = site_mt.aggregate_cols(hl.agg.sum(site_mt.sample_qc.n_snp)) / site_mt.count_cols()
        if i == "MOP":
            print(
                ('Kenya Moi Average Indel and SNP Count:\n\tAverage Indel Count - {:.0f}\n\tAverage SNP count - '
                 '{:.0f}\n').format(
                    (del_avg + ins_avg), snp_avg))
        elif i == "MAP":
            print(
                ('Uganda Average Indel and SNP Count:\n\tAverage Indel Count - {:.0f}\n\tAverage SNP count - {:.0f}'
                 '\n').format(
                    (del_avg + ins_avg), snp_avg))
        elif i == 'KWP':
            print(
                ('Kenya Kemri Average Indel and SNP Count:\n\tAverage Indel Count - {:.0f}\n\tAverage SNP count - '
                 '{:.0f}\n').format(
                    (del_avg + ins_avg), snp_avg))
        elif i == "AAP":
            print(
                ('Ethiopia Average Indel and SNP Count:\n\tAverage Indel Count - {:.0f}\n\tAverage SNP count - '
                 '{:.0f}\n').format(
                    (del_avg + ins_avg), snp_avg))
        elif i == "CTP":
            print(
                ('South Africa Average Indel and SNP Count:\n\tAverage Indel Count - {:.0f}\n\tAverage SNP count - '
                 '{:.0f}\n').format(
                    (del_avg + ins_avg), snp_avg))


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
         '\n').format(
            (del_avg + ins_avg), snp_avg)
    )


# Path for metadata file containing reported sex, collaborator participant id
meta_data = ('gs://neurogap/Pilot_Data_HailQC/sample_AddisEthiopia_KEMRIKenya_MakerereUganda_MoiKenya_UCTSouthAfrica_'
             '.tsv')

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

# Annotating the matrix tables with sample QC data
mt_list = [hl.sample_qc(mt, name='sample_qc') for mt in mt_list]

# Annotating the matrix tables with variant QC data
mt_list = [hl.variant_qc(mt, name='variant_qc') for mt in mt_list]

# Annotating matrix tables with metadata from the meta table (see annotate_meta() for details)
mt_list = [annotateMeta(mt, meta, 'chip_well_barcode') for mt in mt_list]

# Annotating matrix tables with sex filter data (see check_sex() for details)
mt_list = [checkSex(mt) for mt in mt_list]

# Annotating matrix table with sample call rate filter (see check_sample_cr_filter() for details)
mt_list = [checkSampleCrFilter(mt, .98) for mt in mt_list]


# Joining the matrix tables using union_cols(), will need to then to annotate with row data 
for i in range(len(mt_list) - 1):
    if i == 0:
        mt = mt_list[i].union_cols(mt_list[i + 1])
    else:
        mt0 = mt.union_cols(mt_list[i + 1])
        mt = mt0


# Creating a list of site IDs to annotate globals 
siteIDs = hl.array(['MOP', 'MAP', 'KWP', 'AAP ', 'CTP'])
# Annotating globals to the site IDs
mt = mt.annotate_globals(location=siteIDs)


# Creating lists for variant call rate, hwe pval and maf in joined mt
mt = mt.annotate_rows(var_call_rate=hl.empty_array('float64'))
mt = mt.annotate_rows(hwe_pval=hl.empty_array('float64'))
mt = mt.annotate_rows(maf=hl.empty_array('float64'))


# Annotating var_call_rate list with variant call rates
for mt_next in mt_list:
    mt = mt.annotate_rows(var_call_rate=mt.var_call_rate.append(mt_next.index_rows(mt.row_key).variant_qc.call_rate))

# Annotating hwe_pval list with hwe pvals
for mt_next in mt_list:
    mt = mt.annotate_rows(hwe_pval=mt.hwe_pval.append(mt_next.index_rows(mt.row_key).variant_qc.p_value_hwe))

# Annotating maf list with mafs
for mt_next in mt_list:
    mt = mt.annotate_rows(maf=mt.maf.append(hl.min(mt_next.index_rows(mt.row_key).variant_qc.AF)))


# Annotating joined mt with lists for var call rate, hwe pval, and maf flags
mt = mt.annotate_rows(var_cr_flag=hl.empty_array('bool'))
mt = mt.annotate_rows(hwe_pval_flag=hl.empty_array('bool'))
mt = mt.annotate_rows(maf_flag=hl.empty_array('bool'))


# Filling var call rate flag list
for mt_next in mt_list:
    mt = mt.annotate_rows(
        var_cr_flag=mt.var_cr_flag.append(mt_next.index_rows(mt.row_key).variant_qc.call_rate < .95))

# Filling hwe pval flag list
for mt_next in mt_list:
    mt = mt.annotate_rows(
        hwe_pval_flag=mt.hwe_pval_flag.append(mt_next.index_rows(mt.row_key).variant_qc.p_value_hwe <= 1e-03))

# Filling maf flag list
for mt_next in mt_list:
    mt = mt.annotate_rows(maf_flag=mt.maf_flag.append(hl.min(mt_next.index_rows(mt.row_key).variant_qc.AF) <= 0.005))

# Sanity check count
printCount(mt)


# Writing out the matrix table with annotated filter information
mt.write('gs://neurogap/Pilot_Data_HailQC/pilotMT_preQC.mt')

# Reading in the matrix table with all of the site mts combined
mt_joint = hl.read_matrix_table('gs://neurogap/Pilot_Data_HailQC/pilotMT_preQC.mt')


# Filtering the joint dataset into autosomes and X chrom data - will be conducting QC on them separately
mt_auto = mt_joint.filter_rows(mt_joint.locus.contig == "X", keep=False)
mt_x = mt_joint.filter_rows(mt_joint.locus.contig == "X")

# Sanity check counts
printCount(mt_joint)
printCount(mt_auto)
printCount(mt_x)


'''
Getting counts of variants removed by site
'''

mt_rows = mt_joint.rows()

# Fixing variant call rate filter to only one flag, 2 flags was from previous version
mt_rows = mt_rows.annotate(var_cr_flag=mt_rows.var_cr_flag_1)
mt_rows = mt_rows.drop(mt_rows.var_cr_flag_1)
mt_rows = mt_rows.drop(mt_rows.var_cr_flag_2)

var_cr_counts = mt_rows.aggregate(hl.agg.array_agg(lambda x: hl.agg.counter(x), mt_rows.var_cr_flag))

maf_counts = mt_rows.aggregate(hl.agg.array_agg(lambda x: hl.agg.counter(x), mt_rows.maf_flag))

hwe_counts = mt_rows.aggregate(hl.agg.array_agg(lambda x: hl.agg.counter(x), mt_rows.hwe_pval_flag))


# Calculates relatedness using pc_relate for all samples in a matrix table
# Annotates a column which flags those who failed the relatedness filter as True
pc_rel = hl.pc_relate(mt_auto.GT, 0.001, k=10, statistics='kin')
pairs = pc_rel.filter(pc_rel['kin'] > 0.125)
related_samples_to_remove = hl.maximal_independent_set(pairs.i, pairs.j, keep=False)
mt_auto = mt_auto.annotate_cols(related_filter=hl.is_defined(related_samples_to_remove[mt_auto.col_key]))


'''
--- Conducting QC ---

QC Steps:
    snp call rate
    sample call rate
    sex violations
    maf
    hwe
    relatedness
'''

# Getting total counts for samples/variants
printCount(mt_auto)

# Getting variant counts per site
siteVarCount(mt_auto)

# Getting sample counts per site
siteSampleCount(mt_auto)

# Getting the indel counts
hl.summarize_variants(mt_auto)


# SNP call rate 1st pass filtering 
mt_qc = mt_auto.filter_rows((mt_auto.var_cr_flag.contains(True) == True), keep=False)

# Printing out counts post filter
printFilterCounts('SNP call rate', mt_qc, mt_auto.count_rows(), 'variants')

# Getting the indel counts
hl.summarize_variants(mt_qc)

# Getting counts per site for variants and samples post filter
siteVarCount(mt_qc)
siteSampleCount(mt_qc)


# Variable for count pre filter
count = mt_qc.count_cols()

# Sample call rate filtering
mt_qc = mt_qc.filter_cols(mt_qc.sample_cr_filter == True, keep=False)

# Printing out the counts post filter
printFilterCounts('sample call rate', mt_qc, count, 'samples')

# Getting the indel counts
hl.summarize_variants(mt_qc)

# Getting counts per site for variants and samples post filter
siteVarCount(mt_qc)
siteSampleCount(mt_qc)


# This is an implementation of the sex violations part of the check_sex filter adapted from plink
# # Variable for count pre filter
# count = mt_qc.count_cols()

# sexCheck_fail = tmp.cols()
# mt_qc = mt_qc.anti_join_cols(sexCheck_fail)

# # Printing out the counts post filter
# printFilterCounts('sex violations', mt_qc, count, 'samples' )

# # Getting counts per site for variants and samples post filter
# siteVarCount(mt_qc)
# siteSampleCount(mt_qc)

# mt_qc.aggregate_cols(hl.agg.counter(mt_qc.reported_sex))


# Variable for count pre filter
count = mt_qc.count_cols()

# Sex violations filtering
mt_qc = mt_qc.filter_cols(mt_qc.sex_filter == True, keep=False)

# Printing out the counts post filter
printFilterCounts('sex violations', mt_qc, count, 'samples')

# Getting the indel counts
hl.summarize_variants(mt_qc)

# Getting counts per site for variants and samples post filter
siteVarCount(mt_qc)
siteSampleCount(mt_qc)


# Variable for count pre filter
count = mt_qc.count_rows()

# MAF filtering
mt_qc = mt_qc.filter_rows((mt_qc.maf_flag.contains(True) == True), keep=False)

# Printing out counts post filter
printFilterCounts('MAF', mt_qc, count, 'variants')

# Getting the indel counts
hl.summarize_variants(mt_qc)

# Getting counts per site for variants and samples post filter
siteVarCount(mt_qc)
siteSampleCount(mt_qc)


# Variable for count pre filter
count = mt_qc.count_rows()

# HWE filtering
mt_qc = mt_qc.filter_rows((mt_qc.hwe_pval_flag.contains(True) == True), keep=False)

# Printing out counts post filter
printFilterCounts('HWE', mt_qc, count, 'variants')

# Getting the indel counts
hl.summarize_variants(mt_qc)

# Getting counts per site for variants and samples post filter
siteVarCount(mt_qc)
siteSampleCount(mt_qc)


# Variable for count pre filter
count = mt_qc.count_cols()

# Relatedness filtering
mt_qc = mt_qc.filter_cols(mt_qc.related_filter == True, keep=False)

# Printing out the counts post filter
printFilterCounts('relatedness', mt_qc, count, 'samples')

# Checking the indel counts
hl.summarize_variants(mt_qc)

# Getting counts per site for variants and samples post filter
siteVarCount(mt_qc)
siteSampleCount(mt_qc)

# Checking the final counts post QC
printCount(mt_qc)

# Getting counts per site for variants and samples post QC
siteVarCount(mt_qc)
siteSampleCount(mt_qc)


# Writing out the cleaned data as a matrix table
mt_qc.write('gs://neurogap-pilot-clean/NeuroGAP_pilotAuto_clean.mt')

# Writing out the cleaned data as a bgzipped vcf
hl.export_vcf(mt_qc, 'gs://neurogap-pilot-clean/NeuroGAP_pilotAuto_clean.vcf.bgz')


# reading in the post QC dataset to fix/add metrics for output to plink format
mt = hl.read_matrix_table('gs://neurogap/Pilot_Data_HailQC/NeuroGAP_pilotAuto_clean.mt')

# Changing the sex from F and M to True and False so the plink output formats correctly
mt = mt.annotate_cols(sex=hl.if_else((mt.reported_sex == 'F'), True, False))

# Checking that my annotate did as I thought
mt.aggregate_cols(hl.agg.counter(mt.sex == True))

# Reading in table with FID info
fid = hl.import_table('gs://neurogap/Pilot_Data_HailQC/fid_info.csv', delimiter=',')

fid = fid.key_by(fid.siteID)

# Annotating the postQC matrix table with the proper FID information for output to plink
mt = mt.annotate_cols(fid=fid[mt.siteID].FID)

# checking that the annotation worked properly
count = mt.aggregate_cols(hl.agg.counter(mt.fid))

# exporting the final version of the postQC autosomal dataset to plink format
hl.export_plink(mt, 'gs://neurogap-pilot-clean/NeuroGAP_pilotAuto_clean', ind_id=mt.collab_PID, fam_id=mt.fid,
                is_female=mt.sex, varid=mt.rsid)
