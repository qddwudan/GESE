#### This GESE pipeline requires:
python 2.7 (packages re, os, argparse, logging, math)
R (Rscript executable file)
plink2 (plink v1.9 executable file)



### My GESE pipeline includes files:


GESEpipeline.py: starting program, accepting parameter inputs
	input parameters for 'run':
		- project_name: a unique project name for the run, it is used to distinguish different run of the pipeline
                
        - script_dir: the directory for the GESE pipeline

        - Rscript_exec: the location of the Rscript executable file

        - plink_exec: the location of the plink v1.9 executable file

		- plinkfile: the plink file name of the study (You should have the project.bim, project.bed, project.fam file ready for your study data, then put 'project' for this parameter here)

		- studyAnno: the annotated file for the variants in the study using WGSA annotator pipeline.

		- databaseAnno: the annotated file for the variants in the reference data using WGSA annotator pipeline. A good example of this is the ExAC database.

		- fam_file: a file containing the complete pedigree information of the sequenced families. It should include all the subjects in the family even if they were not sequenced. The columns of this file includes at least, FID (family ID), IID (individual ID), faID (father ID, NA if missing), moID(mother ID, NA if missing), and sex.

		- pheno_file: a tab-delimited file containing the phenotypes of the sequenced subjects to be included in the analysis, it should include at least FID, IID, and a column for the disease status in the header. The column name for disease status can be specified using pheno.

		- pheno: specify the column name for disease status in the pheno_file.

		- missense: whether missense variants are considred

		- sift: whether to use sift as a criterion for filtering damaging variants

		- polypheno2: whether to use polyphen2-HVAR as a criterion for filtering damaging variants

        - fathmm: whehter to include detelrious variants predicted by famthmm

		- minCADD: include variants with CADD score greater than this value.

		- maxMAF: include the variants with MAF less than this value, in Europeans in UK10K, 1000GP, and ExAC database.

		- mafMAF_control: include the variants with MAFs less than this value, in the controls in the study 

        - segOnly: whether t ocompute the segregation information only

        - recessive: if segOnly, recessive mode of inheritance is used to compute segregation information

        - dominant: if segOnly, dominant mode of inheritance is used to compute segregation information

        - CH: if segOnly, (pseudo) compound heterozygous mode of inheritance is used to compute segregation information. More information of GESE R package manual.

		- weight_fam: can provide an optional weighting for the families in a file. its dimenstion could be (number of families)x(2). The first column should be family name (column name FID). If the weights for the families are the same for all the genes, the second column should just be weight (columns name "weight"), otherwise the second column and above should be the gene names (columns names are corresponding GENE names). 

	other input parameters for 'rerun' to compute weighted GESE:


		- qpheno: the quantitative phenotype used to weigh the families

        - weight_given: whether the weight matrix has already been created before. If this is false, new weight matrix will be stored in the file with name given by familyWeight parameter.

        - weight_fam: file name of weight matrix, may or may not exist already. 
    
        Note that other parameters required including program_name, script_dir, Rscript_exec, plink_exec, plink_file, studyAnno, databaseAnno, fam_file, pheno_file, pheno should stay the same as running 'run'. Calling 'run' first is required before computing weighted GESE using 'rerun'.'



filter_variants.py: used to filter the variants based on given criterion, made specifically for WGAS Annotator pipeline output

create_subj_list.R: create the subjects that need to be included, this is indicated by the pheno_file (only subjects are included if in this file)

create_variant_list.R: create the variant list that need to be included. It merges the annotation information with the bim file, write the alternative allele file out for recoding to RAW file.

runGESE.R: the file that uses all the input files to run GESE function.

makeWeightMatrix.R: a R function for making the weight matrix using CADD quantitative phenotype information of the subjects.



Example to run the pipeline:
in R
.libPaths("/usr/bin/library/")
install.packages("GESE_1.0.tar.gz", repos=NULL, type="source", INSTALL_opts = c('--no-lock'))

####(LoF or (CADD and SIFT and PPH2)) and MAF filters are included. 
python /udd/redaq/Segregation_analysis/GESE_pipeline/GESEpipeline.py run --project_name ICGN_newF_LoF2_GOLD2 --script_dir /udd/redaq/Segregation_analysis/GESE_pipeline/ --Rscript_exec /udd/redaq/R-3.2.2/bin/Rscript --plink_exec plink2 --plink_file ICGN_EOCOPD_snpsOnly_Final_newFilter --studyAnno /udd/redaq/try/exome1650.annotated.snp --databaseAnno /udd/redaq/try/ExAC.r0.3.annotated.snp --fam_file /udd/redaq/ICGN_EOCOPD/QC/exome1650_4thCleaned.fam --pheno_file ICGN_EOCOPD_extreme_gold234_071916.phe --pheno AFFECT --minCADD -1 --maxMAF 0.001 --maxMAF_control 0.01 > ICGN_newF_LoF2_GOLD2_output.out


### create a matrix function and rerun GESE using the weight (weight depends on phenotype).
### use phenotype only
python /udd/redaq/Segregation_analysis/GESE_pipeline/GESEpipeline.py rerun --project_name ICGN_newF_LoF2_GOLD2 --script_dir /udd/redaq/Segregation_analysis/GESE_pipeline/ --Rscript_exec /udd/redaq/R-3.2.2/bin/Rscript --plink_exec plink2 --plink_file ICGN_EOCOPD_snpsOnly_Final_newFilter --studyAnno /udd/redaq/try/exome1650.annotated.snp --databaseAnno /udd/redaq/try/ExAC.r0.3.annotated.snp --fam_file /udd/redaq/ICGN_EOCOPD/QC/exome1650_4thCleaned.fam --pheno_file ICGN_EOCOPD_extreme_gold234_071916.phe --pheno AFFECT --qpheno residuals --weight_fam ICGN_newF_LoF_GOLD2_weight_phenoOnly.txt > ICGN_newF_LoF2_GOLD2_output.out




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

08/08/2016
- Add in FATHMM as a filtering criteria in filter_variants.py

08/22/2016
- Revised how gene symbol is defined. WGSA0.55 used previous HGNC gene symbol, not consistent with EXAC. If the gene name collumn in WGSA output is not in any of the SNPEFF annotated gene name, then SNPEFF gene name corresponding to the most functionally deleterious annotation is used.
- Add in arguments for the pipeline to specify script location, R and plink version 1.9 executable file locations.

11/15/2016
- In GESE package, updated the condSegProbF function and the segProb function, so each common ancestor is considered separately.
- In GESE package, updated the segProb and the GESE function, gives an MAF estimate for variants absent in the reference database, but present in the study
- Update package description