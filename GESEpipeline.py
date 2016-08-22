import re
import os
import argparse
import logging

from os.path import basename
from os.path import splitext
from filter_variants import select_variants

###### A python script for the pipeline of running GESE

##### Input files
## @ Annotation file for variants in the study obtained using WGSA, WGSA
## @ Annotation file for variants in public database (ExACr0.3)
## @ plink bed, fam, bim file for the study
## @ Complete pedigree info file
## @ Phenotype file

##### Input parameters:
### @ CADD cutoff
### @ SNPEFF categories
### @ MAF in EUR
### @ 



def runGESE(project_name, scriptPATH, RscriptEXE, plinkEXE, plinkfile, studyAnno, databaseAnno, fam_file, pheno_file, pheno, missense=False, sift=False, polyphen2=False, fathmm=False, cadd_min=-1, maf_max=-1, maf_max_control=-1, segOnly = False, recessive = False, dominant = False, CH = False, familyWeight=None):
	
	logger = logging.getLogger(__name__)
    	logger.setLevel(logging.DEBUG)
    	formatter = logging.Formatter("%(levelname)s - %(asctime)-15s %(message)s", "%Y-%m-%d %H:%M:%S")
	ch = logging.StreamHandler()
	ch.setLevel(logging.DEBUG)
	ch.setFormatter(formatter)
	logger.addHandler(ch)

	fh = logging.FileHandler(project_name+'log_pipeline.txt')
	fh.setLevel(logging.DEBUG)
	fh.setFormatter(formatter)
	logger.addHandler(fh)
	
	#### get complete annotation information for databaseAnno
	logger.info("annotation file for the database")
	outputfileAll = "ExAC_r0.3.annotated_anno_all.txt"
	#cmdline = "select_variants "+ databaseAnno + " "+ outputfileAll +" False -1 1 True False False False"
	#logger.info(cmdline)
	#select_variants(databaseAnno, outputfileAll, missense=False, minCADD=-1,  maxMAF=1, noFilter_conseq=True)
	
	
	#### select variants in study Anno
	logger.info("Filter the variants in the study")
	outputfile = splitext(basename(studyAnno))[0] + "_" + project_name + "_anno_filtered.txt"
	#if not os.path.isfile(outputfile):
	cmdline = "select_variants "+studyAnno+" "+outputfile + " "+str(missense)+" " + str(cadd_min)+ " "+str(maf_max)+" "+ str(sift) +" "+ str(polyphen2) + " " + str(fathmm)
	logger.info(cmdline)
	select_variants(studyAnno, outputfile, missense, cadd_min, maf_max, sift=sift, polyphen2=polyphen2, fathmm = fathmm)

	
	#### select variants in database
	logger.info("Filter the variants in the database")
	outputfile2 = splitext(basename(databaseAnno))[0] + "_" + project_name + "_anno_filtered.txt"
	#if not os.path.isfile(outputfile2):
	cmdline = "select_variants "+ databaseAnno + " "+ outputfile2 +" " + str(missense) + " " + str(cadd_min) + " "+ str(maf_max) +" " + str(sift) + " " + str(polyphen2) + " " + str(fathmm)
	logger.info(cmdline)
	select_variants(databaseAnno, outputfile2, missense, cadd_min, maf_max, sift=sift, polyphen2=polyphen2, fathmm = fathmm)
	
	### write the subject list
	subj_output = splitext(basename(pheno_file))[0] + "_" + project_name + "_fidiid.txt"
	control_output = splitext(basename(pheno_file))[0] + "_" + project_name + "_controls.txt"
	subj_list = "{0} {1}create_subj_list.R {2} {3} {4} {5}".format(RscriptEXE, scriptPATH, pheno_file, pheno, subj_output, control_output)
	logger.info(subj_list)
	os.system(subj_list)
	
	#### compute the MAF in the controls 
	freq_output = splitext(basename(pheno_file))[0] + "_" + project_name + "_maf"
	cmd = "{0} --bfile {1} --keep {2} --nonfounders --freq --out {3}".format(plinkEXE, plinkfile, control_output, freq_output)
	logger.info(cmd)
	os.system(cmd)
	
	### write the variant range and create an allele file
	logger.info("Create variant list and allele file")
	idfile = splitext(basename(studyAnno))[0] + "_" + project_name + "_filtered_id.txt"
	allelefile = splitext(basename(studyAnno))[0] + "_" + project_name + "_filtered_allele.txt"
	bimfile = plinkfile+".bim"
	snpInfofile = project_name+"_mapInfo.txt"
	freq_output2 = splitext(basename(pheno_file))[0] + "_" + project_name + "_maf.frq"
	control_removeF = splitext(basename(pheno_file))[0] + "_" + project_name + "_maf_remove.info"
	create_variant_list = "{0} {1}create_variant_list.R {2} {3} {4} {5} {6} {7} {8} {9} {10}".format(RscriptEXE, scriptPATH, outputfile, idfile, bimfile, allelefile, snpInfofile, freq_output2,control_removeF,  maf_max_control, outputfileAll)
	logger.info(create_variant_list)
	os.system(create_variant_list)
	
		
	### extract the variants from the study
	extract_output = project_name +"_filtered"
	extract_output_temp = project_name +"_temp_filtered"
	logger.info("Create raw file for the filtered data")
	## plink2 --bfile eo351_snpFiltId_Final --extract range exome1650_filtered_id.txt --keep eo351_extreme_fidiid.txt --make-bed --out eo351_filtered
	cmd = "{0} --bfile {1} --keep {2} --make-bed --out {3}".format(plinkEXE, plinkfile, subj_output, extract_output_temp)
	#plink2 --bfile eo351_filtered --recode A --recode-allele exome1650_filtered_allele.txt --out eo351_filtered
	logger.info(cmd)
	os.system(cmd)
	cmd = "{0} --bfile {1} --extract {2} --make-bed --out {3}".format(plinkEXE, extract_output_temp, idfile, extract_output)
	logger.info(cmd)
	os.system(cmd)
	
	#### remove the variants with MAF == 0
	poly_output = extract_output + "_poly"
	cmdmaf = "{0} --bfile {1} --nonfounders --reference-allele {2} --maf 0.0000000001 --make-bed --out {3}".format(plinkEXE, extract_output, allelefile, poly_output)
	logger.info(cmdmaf)
	os.system(cmdmaf)
	
	cmd2 = "{0} --bfile {1} --recodeA --recode-allele {2} --out {1}".format(plinkEXE, poly_output, allelefile)
	logger.info(cmd2)
	os.system(cmd2)
	
	### run GESE
	logger.info("Run GESE")
	cmd_gese = "{0} {1}runGESE.R {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12} {13}".format(RscriptEXE, scriptPATH, poly_output, snpInfofile, outputfile2, control_removeF, fam_file, pheno_file, pheno, segOnly, recessive, dominant, CH, familyWeight)
	logger.info(cmd_gese)
	os.system(cmd_gese)
	

def rerunGESE(project_name,  scriptPATH, RscriptEXE, plinkEXE, plinkfile, studyAnno, databaseAnno, fam_file, pheno_file, pheno, qpheno, weightGiven=False, familyWeight=None):
	
	logger = logging.getLogger(__name__)
    	logger.setLevel(logging.DEBUG)
    	formatter = logging.Formatter("%(levelname)s - %(asctime)-15s %(message)s", "%Y-%m-%d %H:%M:%S")
	ch = logging.StreamHandler()
	ch.setLevel(logging.DEBUG)
	ch.setFormatter(formatter)
	logger.addHandler(ch)

	fh = logging.FileHandler(project_name+'log_pipeline.txt')
	fh.setLevel(logging.DEBUG)
	fh.setFormatter(formatter)
	logger.addHandler(fh)
	
	
	#### select variants in study Anno
	#logger.info("Filter the variants in the study")
	outputfile = splitext(basename(studyAnno))[0] + "_" + project_name + "_anno_filtered.tmp"
	#select_variants(studyAnno, outputfile, missense, cadd_min, maf_max, sift=sift, polyphen2=polyphen2)
	
	#### select variants in database
	#logger.info("Filter the variants in the database")
	outputfile2 = splitext(basename(databaseAnno))[0] + "_" + project_name + "_anno_filtered.tmp"
	#select_variants(databaseAnno, outputfile2, missense, cadd_min, maf_max, sift=sift, polyphen2=polyphen2)
	
	### write the subject list
	subj_output = splitext(basename(pheno_file))[0] + "_" + project_name + "_fidiid.txt"
	control_output = splitext(basename(pheno_file))[0] + "_" + project_name + "_controls.txt"
	subj_list = "/udd/redaq/R-3.2.2/bin/Rscript {0}create_subj_list.R {1} {2} {3} {4}".format(scriptPATH, pheno_file, pheno, subj_output, control_output)
	#logger.info(subj_list)
	#os.system(subj_list)
	
	#### compute the MAF in the controls 
	freq_output = splitext(basename(pheno_file))[0] + "_" + project_name + "_maf"
	#cmd = "plink2 --bfile {0} --keep {1} --nonfounders --freq --out {2}".format(plinkfile, control_output, freq_output)
	#logger.info(cmd)
	#os.system(cmd)
	
	### write the variant range and create an allele file
	#logger.info("Create variant list and allele file")
	idfile = splitext(basename(studyAnno))[0] + "_" + project_name + "_filtered_id.txt"
	allelefile = splitext(basename(studyAnno))[0] + "_" + project_name + "_filtered_allele.txt"
	bimfile = plinkfile+".bim"
	snpInfofile = project_name+"_mapInfo.txt"
	freq_output2 = splitext(basename(pheno_file))[0] + "_" + project_name + "_maf.frq"
	control_removeF = splitext(basename(pheno_file))[0] + "_" + project_name + "_maf_remove.info"
	#create_variant_list = "/udd/redaq/R-3.2.2/bin/Rscript {0}create_variant_list.R {1} {2} {3} {4} {5} {6} {7} {8}".format(scriptPATH, outputfile, idfile, bimfile, allelefile, snpInfofile, freq_output2,control_removeF,  maf_max_control)
	#logger.info(create_variant_list)
	#os.system(create_variant_list)
	
		
	### extract the variants from the study
	extract_output = project_name +"_filtered"
	#logger.info("Create raw file for the filtered data")
	#cmd = "plink2 --bfile {0} --extract range {1} --keep {2} --make-bed --out {3}".format(plinkfile, idfile, subj_output, extract_output)
	#logger.info(cmd)
	#os.system(cmd)
	
	#### remove the variants with MAF == 0
	poly_output = extract_output + "_poly"
	#cmdmaf = "plink2 --bfile {0} --nonfounders --reference-allele {1} --maf 0.0000000001 --make-bed --out {2}".format(extract_output, allelefile, poly_output)
	#logger.info(cmdmaf)
	#os.system(cmdmaf)
	
	#cmd2 = "plink --bfile {0} --recodeA --recode-allele {1} --out {0}".format(poly_output, allelefile)
	#logger.info(cmd2)
	#os.system(cmd2)
	
	### create new weights
	weight_file =  familyWeight
	geneSeg_file = poly_output + "_geneSeg.csv"
	varSeg_file =  poly_output + "_varSeg.csv"
	if not weightGiven:
		logger.info("Create weight matrix")
		cmd_r = "/udd/redaq/R-3.2.2/bin/Rscript {0}makeWeightMatrix.R {1} {2} {3} {4} {5} {6} {7} {8} {9}".format(scriptPATH, poly_output, outputfile, outputfile2, geneSeg_file, varSeg_file, pheno_file, pheno, qpheno, weight_file)
		logger.info(cmd_r)
		os.system(cmd_r)
	
	### run GESE
	logger.info("Rerun GESE with weight matrix")
	cmd_gese = "/udd/redaq/R-3.2.2/bin/Rscript {0}runGESE.R {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12}".format(scriptPATH, poly_output, snpInfofile, outputfile2,control_removeF,  fam_file, pheno_file, pheno, False, False, False, False, weight_file)
	logger.info(cmd_gese)
	os.system(cmd_gese)



if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="GESE pipeline")
	subparsers = parser.add_subparsers(help="help")
	run_pipeline = subparsers.add_parser("run", help="run GESE pipeline")
	run_pipeline.set_defaults(which="run")
	run_pipeline.add_argument('--project_name', required=True, help="project name")
	run_pipeline.add_argument('--script_dir', required=True, help="directory for the GESE pipeline files")
	run_pipeline.add_argument('--Rscript_exec', required=True, help="location for the Rscript executable file")
	run_pipeline.add_argument('--plink_exec', required=True, help="location for the plink executable file")
	run_pipeline.add_argument('--plink_file', required=True, help="the plink file (bed, fam, bim files) for the study")
	run_pipeline.add_argument('--studyAnno', required=True, help="the annotation file for the variants in the study")
	run_pipeline.add_argument('--databaseAnno', required=True, help="the annotation file for the variants in the database")
	run_pipeline.add_argument('--fam_file', required=True, help="the file containing complete pedigree info")
	run_pipeline.add_argument('--pheno_file', required=True, help="the file containing the phenotype info")
	run_pipeline.add_argument('--pheno', required=True, help="the column name for the phenotype")
	run_pipeline.add_argument('--missense', action="store_true", help="include missense variants")
	run_pipeline.add_argument('--sift', action="store_true", help="include deleterious variants predicted by sift")
	run_pipeline.add_argument('--polyphen2', action="store_true", help="include deleterious variants predicted by polyphen2")
	run_pipeline.add_argument('--fathmm', action="store_true", help="include deleterious variants predicted by fathmm")
	run_pipeline.add_argument('--minCADD', required=False, help="the minimum CADD value (>)")
	run_pipeline.add_argument('--maxMAF', required=True, help="the maximum MAF (<)")
	run_pipeline.add_argument('--maxMAF_control', required=True, help="the maximum MAF (<) in the unaffected subjects in the data")
	run_pipeline.add_argument('--segOnly', action="store_true", help="Whether to compute the segregation information only")
	run_pipeline.add_argument('--recessive', action="store_true", help="compute segregation information for recessive model")
	run_pipeline.add_argument('--dominant', action="store_true", help="compute segregation information for dominant model")
	run_pipeline.add_argument('--CH', action="store_true", help="compute segregation information for compound heterozygous model")
	
	run_pipeline.add_argument('--weight_fam', required=False, help="a file providing the weighting for the families, could be different for each gene")
	
	run_pipeline2 = subparsers.add_parser("rerun", help="create weight matrix and rerun GESE pipeline")
	run_pipeline2.set_defaults(which="rerun")
	run_pipeline2.add_argument('--project_name', required=True, help="project name")
	run_pipeline2.add_argument('--script_dir', required=True, help="directory for the GESE pipeline files")
	run_pipeline2.add_argument('--Rscript_exec', required=True, help="location for the Rscript executable file")
	run_pipeline2.add_argument('--plink_exec', required=True, help="location for the plink executable file")
	run_pipeline2.add_argument('--plink_file', required=True, help="the plink file for the study")
	run_pipeline2.add_argument('--studyAnno', required=True, help="the annotation file for the variants in the study")
	run_pipeline2.add_argument('--databaseAnno', required=True, help="the annotation file for the variants in the database")
	run_pipeline2.add_argument('--fam_file', required=True, help="the file containing complete pedigree info")
	run_pipeline2.add_argument('--pheno_file', required=True, help="the file containing the phenotype info")
	run_pipeline2.add_argument('--pheno', required=True, help="the column name for the phenotype")
	run_pipeline2.add_argument('--qpheno', required=True, help="the column name for the quantitative phenotype used in the weighting")
	run_pipeline2.add_argument('--weight_given', action="store_true", help="whether the weight is given, create new weights if false")
	run_pipeline2.add_argument('--weight_fam', required=True, help="a file name providing the weighting for the families, could be different for each gene")

	
	args = vars(parser.parse_args())
	if args['which'] == 'run':
		project_name = args['project_name']
		script_dir = args['script_dir']
		Rscript_exec = args['Rscript_exec']
		plink_exec = args['plink_exec']
		plinkfile = args['plink_file']
		studyAnno = args['studyAnno']
		databaseAnno = args['databaseAnno']
		fam_file = args['fam_file']
		pheno_file = args['pheno_file']
		pheno = args['pheno']
		minCADD = args['minCADD']
		maxMAF = args['maxMAF']
		maxMAF_control = args['maxMAF_control']
		familyWeight = args['weight_fam']
		missense = False
		sift = False
		segOnly = False
		polyphen2 = False
		fathmm = False
		recessive = False
		dominant = False
		CH = False
		if parser.parse_args().missense:
			missense = True
		if parser.parse_args().sift:
			sift = True
		if parser.parse_args().polyphen2:
			polyphen2 = True
		if parser.parse_args().segOnly:
			segOnly = True
		if parser.parse_args().dominant:
			dominant = True
		if parser.parse_args().recessive:
			recessive = True
		if parser.parse_args().CH:
			CH = True	
		if parser.parse_args().fathmm:
			fathmm = True
		## set default values for thresholds
		if minCADD == None:
			minCADD = -2
		if maxMAF == None:
			maxMAF = 2
		if maxMAF_control == None:
			maxMAF_control = 2
				
		runGESE(project_name,  script_dir, Rscript_exec, plink_exec, plinkfile, studyAnno, databaseAnno, fam_file, pheno_file, pheno, missense, sift, polyphen2, fathmm, minCADD, maxMAF, maxMAF_control, segOnly, recessive, dominant, CH,  familyWeight)
	elif args['which'] == "rerun":
		project_name = args['project_name']
		plinkfile = args['plink_file']
		studyAnno = args['studyAnno']
		databaseAnno = args['databaseAnno']
		#print anno
		fam_file = args['fam_file']
		pheno_file = args['pheno_file']
		pheno = args['pheno']
		qpheno = args['qpheno']
		familyWeight = args['weight_fam']
		weightGiven = False
		if parser.parse_args().weight_given:
			weightGiven = True
		rerunGESE(project_name,  script_dir, Rscript_exec, plink_exec, plinkfile, studyAnno, databaseAnno, fam_file, pheno_file, pheno, qpheno, weightGiven, familyWeight)

	

	 
