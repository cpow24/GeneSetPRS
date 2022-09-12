# GeneSetPRS
**GENE SET POLYGENIC RISK SCORING** <br />
Exploring the use of enriched gene set polygenic risk scores and machine learning for the prediction of human height <br />

**RUNNING THE CODE** <br />
<br />
**1 -** First, download the individual level data from https://zenodo.org/record/1442755#.YxZpCXbMInI <br />
<br />
**2 -** Then, download the 2018 height GWAS data from https://portals.broadinstitute.org/collaboration/giant/index.php/GIANT_consortium_data_fileBMI and Height GIANT and UK BioBank Meta-analysis Summary Statisticss#BMI_and_Height_GIANT_and_UK_BioBank_Meta-analysis_Summary_Statistics **(Note: the correct files are 'Meta-analysis Wood et al + UKBiobank 2018 GZIP' under 'BMI and Height GIANT and UK BioBank Meta-analysis Summary Statistics'**) - Please rename the GWAS .gz file as 'Height_gwas_2018.txt.gz' <br />
<br />
**3 -** Run the quality control file - 'prs_qc.R' <br />
<br />
**4 -** Run the clumping and thresholding file - 'clump_thresh.R' <br />
<br />
**5 -** Run the LDpred2 and Lassosum file - 'ldpred2_lassosum.R' <br />
<br />
**6 -** Generate gene set scores using 'gene_set_creation.R' <br />
<br />
**7 -** Create necessary gene set files using 'create_gene_set_files.R' <br />
<br />
**8 -** Run the first phase of experiments - 'experiments_1.R' <br />
<br />
**9 -** Run the second stage of experiments - 'experiments_2.ipynb' <br />
<br />
**10 -** Perform the statistical significance testing - 'williams_test.R' <br />
<br />
**11 -** Generate plots - 'plotting.R' <br />
<br />
