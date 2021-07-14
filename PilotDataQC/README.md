This script was created to conduct QC on multiple datasets using the Hail python library and written in a Jupyter Notebook. Data will be read in as a Hail matrix table. See Hail documentation https://hail.is/docs/0.2/methods/impex.html?highlight=import#import for details on what data formats are accepted. This script is written to import vcfs. 

Dependencies: Jupyter Notebook, Python plus Hail 0.2.30 or later

The QC steps and filters used were adapted from Ricopili and Anderson et al. 2010.
There are two separate scripts: One for running QC on the autosomes and PAR region of the X chromosome and one for running QC on the nonPAR region of the X chromosome.
To run both autosomal and X QC the autosomal QC should be run before the X. 

The QC Filtering steps for autosomal and par region of the X chromosome QC are conducted in the following order with these metrics: 
* SNP Call Rate - SNPs with > 5% missingness are removed 
* Sample Call Rate - Individuals with > 2% missingness are removed
* Sex Violations - Individuals whose reported sex do not match genotypic sex are removed
* Minor Allele Frequency - SNPs with MAF < 0.5% are removed
* Hardy Weinberg Equilibrium - SNPs with a HWE p-val < 1e-03 are removed
* Relatedness filtering using PC-Relate with 10 PCs - Individuals with kinship coefficient > .125 removed

The X nonPAR dataset is split into male and female. The male and female datasets then are filtered individually with the following metrics:
* SNP Call Rate - SNPs with > 2% missingness are removed
* Minor Allele Frequency - SNPS with MAF < 1% are removed
* Hardy Weinberg Equilibrium - SNPs with a HWE p-val < 1e-06 are removed

Before removing, variants or samples which fail a QC filter are flagged as True for ease of filtering, and the ability to analyze which QC filters certain datasets are failing.
	
Due to the nature of joining matrix tables in Hail, the variant data for the merged dataset is stored in an array, with n entries, n being the number of original datasets. 
This script is a simplified version of the original, which was used to conduct QC for the 5 NeuroGAP Pilot datasets. 

### Usage: 

If data is stored on google cloud - 
1. Upload this script to your google cloud bucket
2. Start up a cluster using hailctl dataproc start <cluster_name> 
3. Connect to a jupyter notebook using hailctl dataproc connect <cluster_name> notebook
4. Update scripts to work with your data
5. Run Autosomal QC script
6. Run X chromosome QC script   
7. Write out finished QC'd dataset in whatever format is desired.* 

If data is stored on local machine - 
1. Download this script to your local machine 
2. Open a jupyter notebook by typing jupyter notebook into your terminal 
3. Navigate to the script location and open the .ipynb file
4. Update scripts to work with your data
5. Run Autosomal QC script
6. Run X chromosome QC script 
7. Write out finished QC'd dataset in whatever format is desired.*
	
*See https://hail.is/docs/0.2/methods/impex.html?highlight=import#export for more information on what file formats the data can be exported to. 
	
	
### References: 
Anderson, Carl A et al. “Data quality control in genetic case-control association studies.” Nature protocols vol. 5,9 (2010): 1564-73. doi:10.1038/nprot.2010.116
	
Hail Team. Hail 0.2.30-2ae07d872f43. https://github.com/hail-is/hail/releases/tag/0.2.30.
	
Lam, M. et al. RICOPILI: Rapid Imputation for COnsortias PIpeLIne. Bioinformatics https://doi.org/10.1093/bioinformatics/btz633 (2019)
