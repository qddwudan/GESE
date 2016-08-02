#### This GESE pipeline requires:
python 2.7 (packages re, os, argparse, logging, math)
R (Rscript in the PATH)
plink2 (plink v1.9 in PATH)



### My GESE pipeline includes files:


GESEpipeline.py: starting program, accepting parameter inputs
	input parameters for run:
		- project_name: a unique project name for the run, it is used to distinguish different run of the pipeline

		- plinkfile: the plink file name of the study (You should have the project.bim, project.bed, project.fam file ready for your study data, then put 'project' for this parameter here)

		- studyAnno: the annotated file for the variants in the study using WGSA annotator pipeline.

		- databaseAnno: the annotated file for the variants in the reference data using WGSA annotator pipeline. A good example of this is the ExAC database.

		- fam_file: a file containing the complete pedigree information of the sequenced families. It should include all the subjects in the family even if they were not sequenced. The columns of this file includes at least, FID (family ID), IID (individual ID), faID (father ID, NA if missing), moID(mother ID, NA if missing), and sex.

		- pheno_file: a file containing the phenotypes of the sequenced subjects to be included in the analysis, it should include at least FID, IID, and a column for the disease status. The column name for disease status can be specified using pheno.

		- pheno: specify the column name for disease status in the pheno_file.

		- scriptPATH: the directory where the scripts are in

		- missense: whether missense variants are considred

		- sift: whether use sift as a criterion for filtering damaging variants

		- polypheno2: whether use polyphen2-HVAR as a criterion for filtering damaging variants

		- cadd_min: include variants with CADD score greater than this value.

		- maf_max: include the variants with MAF less than this value, in Europeans in UK10K, 1000GP, and ExAC database.

		- maf_max_control: include the variants with MAFs less than this value, in the controls in the study 

		- familyWeight: can provide an optional weighting for the families in a file. its dimenstion could be (number of families)x(2). The first column should be family name (column name FID). If the weights for the families are the same for all the genes, the second column should just be weight (columns name "weight"), otherwise the second column and above should be the gene names (columns names are corresponding GENE names). 

	other input parameters for rerun:

		- qpheno: the quantitative phenotype used to weigh the families


        - weightGiven: whether the weight matrix has already been created before. If this is false, new weight matrix will be stored in the file with name given by familyWeight parameter.

        - familyWeight: file name of weight matrix, may or may not exist already. 



filter_variants.py: used to filter the variants based on given criterion, made specifically for WGAS Annotator pipeline output

create_subj_list.R: create the subjects that need to be included, this is indicated by the pheno_file (only subjects are included if in this file)

create_variant_list.R: create the variant list that need to be included. It merges the annotation information with the bim file, write the alternative allele file out for recoding to RAW file.

runGESE.R: the file that uses all the input files to run GESE function.

makeWeightMatrix.R: a new R function for making the weight matrix using CADD quantitative phenotype information of the subjects.



Example to run the pipeline:
in R
.libPaths("/usr/bin/library/")
install.packages("GESE_1.0.tar.gz", repos=NULL, type="source", INSTALL_opts = c('--no-lock'))

####(LoF or (CADD and SIFT and PPH2)) and MAF filters are included. 
python GESEpipeline.py run --project_name GESE --plink_file eocopd_Final --studyAnno EOCOPD.annotated.snp --databaseAnno ExAC.r0.3.annotated.snp --fam_file pedigree.fam --pheno_file extreme.phe --pheno AFFECT --minCADD 15 --maxMAF 0.001 --sift --polyphen --maxMAF_control 0.01


### create a matrix function and rerun GESE using the weight (weight depends on phenotype).
### use phenotype only
python GESEpipeline.py rerun --project_name rerun_GESE --plink_file eocopd_Final --studyAnno EOCOPD.annotated.snp --databaseAnno try/ExAC.r0.3.annotated.snp --fam_file pedigree.fam --pheno_file extreme.phe --pheno AFFECT --qpheno residuals --weight_fam rerun2_GESE_filtered_poly_weight_phenoOnly.txt



Update log:

05/02/2016
Added maf cutoff in the controls in the data, set default to be 1%. Needs to rerun all the other analysis.

05/03/2016
- Will use the GESE R package
- filter_variants.py is modified, so several possible filtering criterion are included. 
- makeWeightMatrix.R is debugged. Still need to consider what is the best  default CADD value if no variant in the gene is segregating in the family. Right now, we are using min of CADD score for all variants in the gene.
- update trim_oneLineage.R and GESE.R to deal with families with one single subject

05/05/2016
- Update GESE.R getVarSegInfo.R, to include families with one case, but remove singleton families with just one control
- new file trim_unrelated.R, make the function for removing unrelated founder case and families without cases. 


05/09/2016
- Update GESEpipeline.py, create_Variant_list.R, runGESE.R to handle multi-allelic file. 
- also removed the variants with MAF > 0.01 in controls in the data, from the database


05/24/2016
- revised makeWeightMatrix.R, so that families should be trimmed out are removed from the weight matrix. We do not allow missign values in the weight matrix.
- updated GESE R package. Allow different number of families in the weight matrix and our raw file.
- make sure the weight column has name "weight" in getPvalue_resampling.R

06/06/2016
- in GESE package, update GESE.R such that the calculation of segregation information is more efficient, the gene-based and variant based segregation information are also returned.
- remove the two functions getGeneSegInfo and getVarSegInfo. All segregation information can be obtained using GESE function specifying onlySeg=TRUE.
- in runGESE.R, we do not need to compute variant-based and gene-based segregation information separately.
- make sure the weight column has name "weight" in getPvalue_resampling.R

06/23/2016
- In the GESE package, we have very strict checking on segregation. If there is missing genotype in either cases or controls, there would be no segregation event. This behavior may be modified later.

06/24/2016
- In GESE package, added an imputation step, where all missing genotypes were imputed to be homozygous major allele genotypes (0).
- Added getSegInfo.R for computing segregation information for different mode of inheritance.
- In GESE pipeline, added in runGESE.R to get segregation information for other mode of inheritance.

07/19/2016
- Add in an option in GESEpipeline.py, so user can request to compute segregation information only, no GESE test was done.

07/25/2016
- Remove genes with no coverage in the database from results, probably due to coverage difference

