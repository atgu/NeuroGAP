#!/usr/bin/env python
# coding: utf-8

# In[1]:

import hail as hl
hl.init()

# In[2]:

from hail.plot import show
from pprint import pprint
hl.plot.output_notebook()

# In[3]:

'''
Annotating Functions
'''
# Conducts sex imputation statistics for a site. Returns an mt with an annotated 
# imputed sex column & a column which flags those who failed sex filter as True
def checkSex(mt):
    new_mt = hl.impute_sex(mt.GT)
    mt = mt.annotate_cols(imputedSex = new_mt[mt.s])
    return(mt.annotate_cols(sex_filter = mt.imputedSex.is_female != (mt.reported_sex=='F')))

# Annotates the mt with the results of the sample call rate filter
# Creates a column which flags those who failed the filter as True
def checkSampleCrFilter(mt, sample_geno_filter):
    return(mt.annotate_cols(sample_cr_filter = mt.sample_qc.call_rate < sample_geno_filter))

# Annotates a matrix table with the following meta data: collaborator participant id, site id, and 
# autocall call rate. Assumes that meta and mt don't have the same key, keys meta by chip_well_barcode
def annotateMeta(mt,meta,key_by):
    
    # Setting the key of metaDat as chip_well_barcode
    meta = meta.key_by(key_by)
    
    # Adding the reported sex column from the sample data table to the matrix table
    mt = mt.annotate_cols(reported_sex = meta[mt.s].reported_gender)

    # Annotating the matrix table with the projID from meta data
    mt = mt.annotate_cols(collab_PID = meta[mt.s].collaborator_participant_id)

    # Creating a new column in mt for siteID which is the 3 letters from the collab participant id
    mt = mt.annotate_cols(siteID = mt.col.collab_PID[0:3])
     
    # Switching lowercase letters to capital in collabPID and siteID
    mt = mt.annotate_cols(collab_PID = mt.collab_PID.upper())
    mt = mt.annotate_cols(siteID = mt.siteID.upper())

    # Filtering any NA12878 individuals from the matrix table
    samples_to_remove = {'NA12878'}
    set_to_remove = hl.literal(samples_to_remove)
    return(mt.filter_cols(~set_to_remove.contains(mt['collab_PID'])))

# In[4]:

'''
Printing Functions
'''
# Prints out total sample/variant count for an mt
def printCount(mt):
    n = mt.count()
    print('Number of Samples: {}\nNumber of Variants: {}'.format(n[1],n[0]))    

# Prints how many remaining samples/variants there are after filtering. 
# Where mt is the mt post filtering, count is mt.count() pre filter and filter_type is 'samples' or 'variants'
def printFilterCounts(filter_name, mt, count, filter_type):
    if filter_type == 'variants':
        print('After %s filter, %d/%d %s remain.\n' % (filter_name,mt.count_rows(),count,filter_type))
    elif filter_type == 'samples':
        print('After %s filter, %d/%d %s remain.\n' % (filter_name,mt.count_cols(),count,filter_type))

# In[5]:

# Path for metadata file containing reported sex, collaborator participant id
# This file is necessary if original datasets do not contain reported sex info, which is needed to sex violations filtering
meta_data = 'file/path'

# In[6]:
# MT paths is a list of file paths for each of the datasets to be merged and QC'd  
mt_paths = ['file/path1','file/path2']

# In[7]:

# Reading in and creating a list of all of the site matrix tables 
mt_list = [hl.import_vcf(mt_path,force_bgz = True) for mt_path in mt_paths]

# Importing the metadata file as a hail table
meta = hl.import_table(meta_data)

# In[8]:

# Annotating the matrix tables with sample QC data
mt_list = [hl.sample_qc(mt, name = 'sample_qc') for mt in mt_list]

# In[89]:

# Annotating the matrix tables with variant QC data
mt_list = [hl.variant_qc(mt, name = 'variant_qc') for mt in mt_list]

# In[90]:

# Annotating matrix tables with metadata from the meta table (see annotateMeta for details) 
mt_list = [annotateMeta(mt,meta,'chip_well_barcode') for mt in mt_list]

# In[91]:

# Annotating matrix tables with sex filter results (see checkSex for details)
mt_list = [checkSex(mt) for mt in mt_list]

# In[92]:

# Annotating matrix table with sample call rate filter results (see checkSampleCrFilter for details)
mt_list = [checkSampleCrFilter(mt,.98) for mt in mt_list]

# In[93]:

# Joining the matrix tables using union_cols(), will need to then to annotate row arrays with joint row data 
for i in range(len(mt_list)-1):
    if i == 0:
         mt = mt_list[i].union_cols(mt_list[i+1])
    else: 
        mt0 = mt.union_cols(mt_list[i+1])
        mt = mt0

# In[94]:

# Creating a list of site IDs to annotate globals - for the purpose of keeping track of the order of variant data in variant arrays 
siteIDs = hl.array(['ID', 'ID', 'ID'])

mt = mt.annotate_globals(location = siteIDs) 

# In[95]:

# Creating lists for variant call rate, hwe pval and maf in joined mt
mt = mt.annotate_rows(var_call_rate = hl.empty_array('float64'))
mt = mt.annotate_rows(hwe_pval = hl.empty_array('float64'))
mt = mt.annotate_rows(maf = hl.empty_array('float64'))

# In[96]:

# Annotating var_call_rate list with variant call rates for each location
for mt_next in mt_list: 
    mt = mt.annotate_rows(var_call_rate = mt.var_call_rate.append(mt_next.index_rows(mt.row_key).variant_qc.call_rate))

# In[98]:

# Annotating hwe_pval list with hwe pvals for each location
for mt_next in mt_list:
    mt = mt.annotate_rows(hwe_pval = mt.hwe_pval.append((mt_next.index_rows(mt.row_key).variant_qc.p_value_hwe)))

# In[99]:

# Annotating maf list with mafs for each location
for mt_next in mt_list:
    mt = mt.annotate_rows(maf = mt.maf.append(hl.min(mt_next.index_rows(mt.row_key).variant_qc.AF)))

# In[100]:

# Annotating joined mt with lists for var call rate, hwe pval and maf flags
# Each list will contain the respective metric for each of the 5 locations, in the order specified by the global locations array
mt = mt.annotate_rows(var_cr_flag = hl.empty_array('bool'))
mt = mt.annotate_rows(hwe_pval_flag = hl.empty_array('bool'))
mt = mt.annotate_rows(maf_flag = hl.empty_array('bool'))

# In[101]:

# Filling var call rate flag list
for mt_next in mt_list:
    mt = mt.annotate_rows(var_cr_flag = mt.var_cr_flag.append((mt_next.index_rows(mt.row_key).variant_qc.call_rate)<.95))

# In[102]:

# Filling hwe pval flag list
for mt_next in mt_list:
    mt = mt.annotate_rows(hwe_pval_flag = mt.hwe_pval_flag.append((mt_next.index_rows(mt.row_key).variant_qc.p_value_hwe)<=1e-03))

# In[103]:

# Filling maf flag list
for mt_next in mt_list:
    mt = mt.annotate_rows(maf_flag = mt.maf_flag.append(hl.min(mt_next.index_rows(mt.row_key).variant_qc.AF) <= 0.005))

# In[25]:

# Writing out the matrix table with annotated filter information
# Writing out the joint matrix makes the following QC steps run faster
mt.write('path/for/joint/matrix')

# In[9]:

# Reading in the matrix table with all of the site mts combined
joint_data = 'path/for/joint/matrix'
mt_joint = hl.read_matrix_table(joint_data)

# In[12]:

# Calculates relatedness using pc_relate for all samples in a matrix table
# Annotates a column which flags those who failed the relatedness filter as True

pc_rel = hl.pc_relate(mt_joint.GT, 0.001, k=10, statistics='kin')
pairs = pc_rel.filter(pc_rel['kin'] > 0.125)
related_samples_to_remove = hl.maximal_independent_set(pairs.i, pairs.j,keep=False)
mt_joint = mt_joint.annotate_cols(related_filter = hl.is_defined(related_samples_to_remove[mt_joint.col_key]))

# In[11]:

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

# In[128]:

# SNP call rate filtering 
mt_qc = mt_joint.filter_rows((mt_joint.var_cr_flag.contains(True) == True), keep = False)

# Printing out counts post filter
printFilterCounts('SNP call rate', mt_qc, mt_joint.count_rows(), 'variants' )

# In[60]:

# Variable for count pre filter
count = mt_qc.count_cols()

# Sample call rate filtering
mt_qc = mt_qc.filter_cols(mt_qc.sample_cr_filter == True, keep=False)

# Printing out the counts post filter
printFilterCounts('sample call rate', mt_qc, count, 'samples' )

# In[61]:

# Variable for count pre filter
count = mt_qc.count_cols()

# Sex violations filtering
mt_qc = mt_qc.filter_cols(mt_qc.sex_filter == True, keep=False)

# Printing out the counts post filter
printFilterCounts('sex violations', mt_qc, count, 'samples' )

# In[63]:

# Variable for count pre filter
count = mt_qc.count_rows() 

# MAF filtering
mt_qc = mt_qc.filter_rows((mt_qc.maf_flag.contains(True) == True), keep = False)

# Printing out counts post filter
printFilterCounts('MAF', mt_qc, count, 'variants')

# In[64]:

# Variable for count pre filter
count = mt_qc.count_rows() 

# HWE filtering
mt_qc = mt_qc.filter_rows((mt_qc.hwe_pval_flag.contains(True) == True), keep = False)

# Printing out counts post filter
printFilterCounts('HWE', mt_qc, count, 'variants' )

# In[65]:

# Variable for count pre filter
count = mt_qc.count_cols()

# Relatedness filtering
mt_qc = mt_qc.filter_cols(mt_qc.related_filter == True, keep=False)

# Printing out the counts post filter
printFilterCounts('relatedness', mt_qc, count, 'samples')
