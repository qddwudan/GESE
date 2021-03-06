import re
import sys
import os
import math
from os.path import basename
from os.path import splitext
###### A python script for filtering variants on annotation file generated using WGSA

##### Input files
## @ Annotation file for variants obtained using WGSA
## @ output file name
## @ minimum CADD score (>)
## @ maximum MAF (<)
## @ SIFT prediction is D (score < 0.05)
## @ Polypheno2 HVAR prediction is D


def select_variants(annoFile, output, missense, minCADD,maxMAF, noFilter_conseq=False, sift=False, polyphen2=False, fathmm=False):
	
	outlong = splitext(basename(output))[0]+ ".tmp"
	out = open(output, "w")
	out2 = open(outlong, "w")
	out.write("CHR\tPOS\tREF\tALT\tMAF_EXAC\tGENE\n")
	out2.write("CHR\tPOS\tREF\tALT\tMAF_1000G\tMAF_EXAC\tMAF_UK10K\tGENE\tCADD\tCONSEQUENCE\tSIFT\tPOLYPHEN2_HVAR\tFATHMM_score\tFATHMM_PRED\n")
	highmodEffect = ['frameshift_variant', 'splice_acceptor_variant', 'splice_donor_variant', 'stop_lost', 'start_lost', 'stop_gained', 'missense_variant', 'rare_amino_acid_variant', "chromosome", "exon_loss_variant"]
	highEffect = ['frameshift_variant', 'splice_acceptor_variant', 'splice_donor_variant', 'stop_lost', 'start_lost', 'stop_gained', 'rare_amino_acid_variant', "chromosome", "exon_loss_variant"]
	modEffect = ["coding_sequence_variant", "inframe_insertion", "disruptive_inframe_insertion", "inframe_deletion", "disruptive_inframe_deletion", "missense_variant", "splice_region_variant" ]
	#modifierEffect = ["coding_sequence_variant", "downstream_gene_variant", "exon_variant", "gene_variant", "intergenic_region", "conserved_intergenic_variant", "intragenic_variant", "intron_variant", "conserved_intron_variant", "miRNA", "transcript_variant", "regulatory_region_variant", "upstream_gene_variant", "3_prime_UTR_variant", "5_prime_UTR_variant"]
	otherEffect = ['synonymous_variant', "coding_sequence_variant","exon_variant", "transcript_variant",'initiator_codon_variant', 'stop_retained_variant', 'start_retained', 'stop_retained_variant']
	
	minCADD = float(minCADD)
	maxMAF = float(maxMAF)
	with open(annoFile, 'r') as f:
		for line in f:
			
			#chr=pos=ref=alt=cadd_index=maf_1000G_index=maf_ExAC_NFE_index=maf_UK10K_index=snpeff_index=gene_index=0
			if "#" in line:
				header = line.split("\t")
				chr_index = header.index("#chr")
				pos_index = header.index("pos")
				ref_index = header.index("ref")
				alt_index = header.index("alt")	
				
				cadd_index = header.index("CADD_phred")
				fathmm_index = header.index("FATHMM_score")

				maf_1000G_index = header.index("1000Gp3_EUR_AF")
				maf_ExAC_NFE_index = header.index("ExAC_NFE_AF")
				maf_UK10K_index = header.index("TWINSUK_AF")
				snpeff_index = header.index("SnpEff_refseq_summary")
				gene_index = header.index("genename")
				sift_index = header.index("SIFT_pred")
				pphen2_index = header.index("Polyphen2_HVAR_pred")
				fathmm_pred_index = header.index("FATHMM_pred")
				continue
			
			fields = line.split("\t")
			chr = fields[chr_index]
			pos = fields[pos_index]
			ref = fields[ref_index]
			alt = fields[alt_index]
			cadd = fields[cadd_index]
			fathmmD = fields[fathmm_index]
			maf_1000G = fields[maf_1000G_index]
			maf_ExAC_NFE = fields[maf_ExAC_NFE_index]
			maf_UK10K = fields[maf_UK10K_index]
			snpeff = fields[snpeff_index]
			siftD = fields[sift_index]
			pphen2 = fields[pphen2_index]
			fathmm_pred = fields[fathmm_pred_index]
			#### find gene name
			annoInfo = snpeff.split("|")
			geneSet = [x.split("(")[0] for x in annoInfo]
			gene = fields[gene_index]  #annoInfo[0].split("(")[0]
			
			highInfo = [any([y in x for y in highEffect]) for x in annoInfo]
			highmodInfo = [any([y in x for y in highmodEffect]) for x in annoInfo]
			modInfo = [any([y in x for y in modEffect]) for x in annoInfo]
			otherInfo = [any([y in x for y in otherEffect]) for x in annoInfo]
			
			if "|" in gene or ";" in gene or gene == "." or gene not in geneSet:
				infoTemp = annoInfo[0]					
				if any(otherInfo):
					infoTemp = [annoInfo[x] for x in range(len(otherInfo)) if otherInfo[x]][0]
				if any(modInfo):
					infoTemp = [annoInfo[x] for x in range(len(modInfo)) if modInfo[x]][0]
				if any(highmodInfo):
					infoTemp = [annoInfo[x] for x in range(len(highmodInfo)) if highmodInfo[x]][0]
				if any(highInfo):
					infoTemp = [annoInfo[x] for x in range(len(highInfo)) if highInfo[x]][0]
				if "(" in infoTemp:
					gene = infoTemp.split("(")[0]
			
			
				
						
			if maf_1000G == '.':
				maf_1000G = 'nan'
			if maf_ExAC_NFE == '.':
				maf_ExAC_NFE = 'nan'
			if maf_UK10K == '.':
				maf_UK10K = 'nan'
			
			passFATHMM = True	
			### fathmm is only a filter when the tag is present
			if fathmm:
				if 'D' in fathmm_pred:
					passFATHMM = True
				else:
					passFATHMM = False
			
			
					
				
			if noFilter_conseq:
				login = 'no filter is used \n'
				conseq = snpeff
				out2.write(str(chr) + '\t' + str(pos) + '\t' + str(ref) + '\t' + str(alt) + '\t' + str(maf_1000G) + '\t' + str(maf_ExAC_NFE)	+ '\t' + str(maf_UK10K) + '\t' + str(gene)+ '\t' + str(cadd)+ '\t' + str(conseq) + '\t' + str(siftD) + '\t' + str(pphen2) + '\t' + str(fathmmD) + '\t' + str(fathmm_pred) + '\n')			
				out.write(str(chr)  + '\t' + str(pos) + '\t' + str(ref) + '\t' + str(alt) + '\t' + str(maf_ExAC_NFE) + '\t' + str(gene) + '\n')
			elif not missense and not sift and not polyphen2: ## only LoF  and CADD and FATHMM and MAF threshold included
				login = 'only LoF variants, and CADD, FATHMM and MAF filters are included. \n'
				if float(cadd) > minCADD  and passFATHMM and (math.isnan(float(maf_1000G)) or float(maf_1000G) < maxMAF)  and ( math.isnan(float(maf_UK10K)) or float(maf_UK10K) < maxMAF) and (math.isnan(float(maf_ExAC_NFE)) or float(maf_ExAC_NFE) < maxMAF) and any(highInfo):
					conseq = '|'.join([x for x in highEffect if x in snpeff])
					out2.write(str(chr) + '\t' + str(pos) + '\t' + str(ref) + '\t' + str(alt) + '\t' + str(maf_1000G) + '\t' + str(maf_ExAC_NFE)	+ '\t' + str(maf_UK10K) + '\t' + str(gene)+ '\t' + str(cadd)+ '\t' + str(conseq)  + '\t' + str(siftD) + '\t' + str(pphen2) + '\t' + str(fathmmD) + '\t' + str(fathmm_pred) +'\n')			
					out.write(str(chr)  + '\t' + str(pos) + '\t' + str(ref) + '\t' + str(alt) + '\t' + str(maf_ExAC_NFE) + '\t' + str(gene) + '\n')
			elif missense and not sift and not polyphen2: ## only LoF (or also missense) and CADD and FATHMM and MAF included
				login = 'only LoF and missense variants, and CADD,FATHMM and MAF filters are included. \n'
				if float(cadd) > minCADD  and passFATHMM and (math.isnan(float(maf_1000G)) or float(maf_1000G) < maxMAF)  and ( math.isnan(float(maf_UK10K)) or float(maf_UK10K) < maxMAF) and (math.isnan(float(maf_ExAC_NFE)) or float(maf_ExAC_NFE) < maxMAF) and any(highmodInfo):
					conseq = '|'.join([x for x in highmodEffect if x in snpeff])
					out2.write(str(chr) + '\t' + str(pos) + '\t' + str(ref) + '\t' + str(alt) + '\t' + str(maf_1000G) + '\t' + str(maf_ExAC_NFE)	+ '\t' + str(maf_UK10K) + '\t' + str(gene)+ '\t' + str(cadd)+ '\t' + str(conseq)  + '\t' + str(siftD) + '\t' + str(pphen2) + '\t' + str(fathmmD) + '\t' + str(fathmm_pred) +'\n')			
					out.write(str(chr)  + '\t' + str(pos) + '\t' + str(ref) + '\t' + str(alt) + '\t' + str(maf_ExAC_NFE) + '\t' + str(gene) + '\n')
			elif not missense and sift and polyphen2: ## (LoF or (CADD and SIFT and PPH2)) and MAF included
				login = 'LoF variants, or (any variants that pass CADD and FATHMM threshold and predicted to be damaging by SIFT and PPH2, and MAF filters are included. \n'
				if ((math.isnan(float(maf_1000G)) or float(maf_1000G) < maxMAF)  and ( math.isnan(float(maf_UK10K)) or float(maf_UK10K) < maxMAF) and (math.isnan(float(maf_ExAC_NFE)) or float(maf_ExAC_NFE) < maxMAF)) and (any(highInfo) or (float(cadd) > minCADD  and passFATHMM and ('D' in siftD) and ('D' in pphen2))):
					conseq = '|'.join([x for x in highEffect if x in snpeff])
					out2.write(str(chr) + '\t' + str(pos) + '\t' + str(ref) + '\t' + str(alt) + '\t' + str(maf_1000G) + '\t' + str(maf_ExAC_NFE)	+ '\t' + str(maf_UK10K) + '\t' + str(gene)+ '\t' + str(cadd)+ '\t' + str(conseq) + '\t' + str(siftD) + '\t' + str(pphen2) + '\t' + str(fathmmD) + '\t' + str(fathmm_pred) + '\n')			
					out.write(str(chr)  + '\t' + str(pos) + '\t' + str(ref) + '\t' + str(alt) + '\t' + str(maf_ExAC_NFE) + '\t' + str(gene) + '\n')
			elif missense and sift and polyphen2: ## (LoF or (CADD and SIFT and PPH2 and missense)) and MAF included
				login = 'LoF variants , or missense variants that pass CADD and FATHMM threshold and predicted to be damaging by SIFT and PPH2, and MAF filters are included. \n'
				if ((math.isnan(float(maf_1000G)) or float(maf_1000G) < maxMAF)  and ( math.isnan(float(maf_UK10K)) or float(maf_UK10K) < maxMAF) and (math.isnan(float(maf_ExAC_NFE)) or float(maf_ExAC_NFE) < maxMAF)) and (any(highInfo) or (float(cadd) > minCADD  and passFATHMM and ('D' in siftD) and ('D' in pphen2) and any(highmodEffect))):
					conseq = '|'.join([x for x in highmodEffect if x in snpeff])
					out2.write(str(chr) + '\t' + str(pos) + '\t' + str(ref) + '\t' + str(alt) + '\t' + str(maf_1000G) + '\t' + str(maf_ExAC_NFE)	+ '\t' + str(maf_UK10K) + '\t' + str(gene)+ '\t' + str(cadd)+ '\t' + str(conseq)  + '\t' + str(siftD) + '\t' + str(pphen2) + '\t' + str(fathmmD) + '\t' + str(fathmm_pred) +'\n')			
					out.write(str(chr)  + '\t' + str(pos) + '\t' + str(ref) + '\t' + str(alt) + '\t' + str(maf_ExAC_NFE) + '\t' + str(gene) + '\n')
			elif not missense and sift and not polyphen2: ## (LoF or (CADD and SIFT)) and MAF included
				login = '(LoF variants or (any variants that pass CADD and FATHMM thresholds, and predicted to be damaging by SIFT, and MAF filters are included. \n'
				if ((math.isnan(float(maf_1000G)) or float(maf_1000G) < maxMAF)  and ( math.isnan(float(maf_UK10K)) or float(maf_UK10K) < maxMAF) and (math.isnan(float(maf_ExAC_NFE)) or float(maf_ExAC_NFE) < maxMAF)) and (any(highInfo) or (float(cadd) > minCADD and passFATHMM and ('D' in siftD))):
					conseq = '|'.join([x for x in highEffect if x in snpeff])
					out2.write(str(chr) + '\t' + str(pos) + '\t' + str(ref) + '\t' + str(alt) + '\t' + str(maf_1000G) + '\t' + str(maf_ExAC_NFE)	+ '\t' + str(maf_UK10K) + '\t' + str(gene)+ '\t' + str(cadd)+ '\t' + str(conseq)  + '\t' + str(siftD) + '\t' + str(pphen2) + '\t' + str(fathmmD) + '\t' + str(fathmm_pred) +'\n')			
					out.write(str(chr)  + '\t' + str(pos) + '\t' + str(ref) + '\t' + str(alt) + '\t' + str(maf_ExAC_NFE) + '\t' + str(gene) + '\n')
			elif not missense and not sift and polyphen2: ## (LoF or (CADD and polyphen2)) and MAF included
				login = '(LoF variants or (any variants that pass CADD and FATHMM thresholds, and predicted to be damaging by Polyphen2, and MAF filters are included. \n'
				if ((math.isnan(float(maf_1000G)) or float(maf_1000G) < maxMAF)  and ( math.isnan(float(maf_UK10K)) or float(maf_UK10K) < maxMAF) and (math.isnan(float(maf_ExAC_NFE)) or float(maf_ExAC_NFE) < maxMAF)) and (any(highInfo) or (float(cadd) > minCADD and passFATHMM and ('D' in pphen2))):
					conseq = '|'.join([x for x in highEffect if x in snpeff])
					out2.write(str(chr) + '\t' + str(pos) + '\t' + str(ref) + '\t' + str(alt) + '\t' + str(maf_1000G) + '\t' + str(maf_ExAC_NFE)	+ '\t' + str(maf_UK10K) + '\t' + str(gene)+ '\t' + str(cadd)+ '\t' + str(conseq)  + '\t' + str(siftD) + '\t' + str(pphen2) + '\t' + str(fathmmD) + '\t' + str(fathmm_pred) +'\n')			
					out.write(str(chr)  + '\t' + str(pos) + '\t' + str(ref) + '\t' + str(alt) + '\t' + str(maf_ExAC_NFE) + '\t' + str(gene) + '\n')
			elif missense and sift and not polyphen2: ## (LoF or (CADD and SIFT and missense)) and MAF included
				login = 'LoF variants, or missense variants that pass CADD and FATHMM thresholds, and predicted to be damaging by SIFT, and MAF filters are included. \n'
				if ((math.isnan(float(maf_1000G)) or float(maf_1000G) < maxMAF)  and ( math.isnan(float(maf_UK10K)) or float(maf_UK10K) < maxMAF) and (math.isnan(float(maf_ExAC_NFE)) or float(maf_ExAC_NFE) < maxMAF)) and (any(highInfo) or (float(cadd) > minCADD and passFATHMM and ('D' in siftD) and any(highmodEffect))):
					conseq = '|'.join([x for x in highmodEffect if x in snpeff])
					out2.write(str(chr) + '\t' + str(pos) + '\t' + str(ref) + '\t' + str(alt) + '\t' + str(maf_1000G) + '\t' + str(maf_ExAC_NFE)	+ '\t' + str(maf_UK10K) + '\t' + str(gene)+ '\t' + str(cadd)+ '\t' + str(conseq)  + '\t' + str(siftD) + '\t' + str(pphen2) + '\t' + str(fathmmD) + '\t' + str(fathmm_pred) +'\n')			
					out.write(str(chr)  + '\t' + str(pos) + '\t' + str(ref) + '\t' + str(alt) + '\t' + str(maf_ExAC_NFE) + '\t' + str(gene) + '\n')
			elif missense and not sift and polyphen2: ## (LoF or (CADD and polyphen2 and missense)) and MAF included
				login = 'LoF variants, or missense variants that pass CADD and FATHMM thresholds, and predicted to be damaging by Polyphen2, and MAF filters are included.  \n'
				if ((math.isnan(float(maf_1000G)) or float(maf_1000G) < maxMAF)  and ( math.isnan(float(maf_UK10K)) or float(maf_UK10K) < maxMAF) and (math.isnan(float(maf_ExAC_NFE)) or float(maf_ExAC_NFE) < maxMAF)) and (any(highInfo) or (float(cadd) > minCADD and passFATHMM  and ('D' in pphen2) and any(highmodEffect))):
					conseq = '|'.join([x for x in highmodEffect if x in snpeff])
					out2.write(str(chr) + '\t' + str(pos) + '\t' + str(ref) + '\t' + str(alt) + '\t' + str(maf_1000G) + '\t' + str(maf_ExAC_NFE)	+ '\t' + str(maf_UK10K) + '\t' + str(gene)+ '\t' + str(cadd)+ '\t' + str(conseq)  + '\t' + str(siftD) + '\t' + str(pphen2) + '\t' + str(fathmmD) + '\t' + str(fathmm_pred) +'\n')			
					out.write(str(chr)  + '\t' + str(pos) + '\t' + str(ref) + '\t' + str(alt) + '\t' + str(maf_ExAC_NFE) + '\t' + str(gene) + '\n')

			
			
	print login
	out.close()
	out2.close()

if __name__ == "__main__":
	#udd/redaq/try/exome1650.annotated.snp exome1650.annotated_ICGN_newF_segOnly_01_noCSP2_anno_filtered.txt True -1 0.1 False False False
	annoFile = "/udd/redaq/try/exome1650.annotated.snp"
	output = "exome1650.annotated_ICGN_newF_segOnly_01_noCSP2_anno_filtered.txt"
	minCADD = -1
	maxMAF = 0.1
	missense = True
	sift = False
	polyphen2 = False
	fathmm = False
	#select_variants(annoFile, output, missense, minCADD, maxMAF, noFilter_conseq, sift, polyphen2)
	missense = sys.argv[3] == "True"
	noFilter_conseq = sys.argv[6] == "True"
	sift = sys.argv[7] == "True"
	polyphen2 = sys.argv[8] == "True"
	fathmm = sys.argv[9] == "True"
	select_variants(sys.argv[1], sys.argv[2], missense, sys.argv[4],  sys.argv[5],  noFilter_conseq, sift, polyphen2, fathmm  )
