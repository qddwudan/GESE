##### This file is used to call the main GESE function

#### input
## @ extract_output: the plink raw file name (without .raw), with the bed, bim, fam files
## @ studyAnno: the studyAnno file
## @ databaseAnno: the database annotation file, created in previous steps.
## @ variant_remove: the variants removed due to MAF in controls. It should contain a column chrPos, which is in the format of chr:position. Eg. 1:24523410. This file was also created in the pipeline.
## @ fam_file: the file containing complete pedigree info. It should contain columns FID, IID, faID, moID, sex.
## @ pheno_file: the phenotype information, it should contain the FID, IID and phenotype columns
## @ familyWeightF: a data frame that contains the weighting of the families. The number of rows should equal to the number of families. There should be two columns: FID and weight.

args <- commandArgs(TRUE)
extract_output = args[1]
studyAnno = args[2]
databaseAnno = args[3]
variant_remove = args[4]
fam_file = args[5]
pheno_file  = args[6]
pheno = args[7]
segOnly = args[8]
recessive = args[9]
dominant = args[10]
CH = args[11]
familyWeightF = args[12]


options(stringsAsFactors=FALSE)

### pedigree info
pednew <- read.table(fam_file, header=T)
#set library path .libPaths("")
.libPaths("/udd/redaq/bin/R-2.15.0/library/")
library(kinship2)
library(GESE)
ped <- pedigree(id=pednew$IID, dadid=pednew$faID, momid = pednew$moID, sex=pednew$sex, famid=pednew$FID)

### map info
map <- read.table(paste(extract_output, ".bim", sep=""))
snpInfo <- read.table(studyAnno, header=T)
mapInfo <- merge(map, snpInfo, by.x="V2", by.y="SNP", sort=F)
mapInfo <- mapInfo[,c("V2", "GENE")]
colnames(mapInfo) <- c("SNP", "GENE")

## database info
database = read.table(databaseAnno, sep="\t", header=T)
database$chrPos = paste(database$CHR, database$POS, sep=":")
### have to collapse MAF for multi-allelic variants
mafDatabase = aggregate(x=database$MAF_EXAC, by=list(database$chrPos), FUN="sum", simplify=TRUE, na.rm=T)
db1 = database[!duplicated(database$chrPos),]
database2 = merge( mafDatabase, db1, by.y="chrPos", by.x="Group.1", sort=F, all.x=T)
varRemove = read.table(variant_remove, header=T)
database3 = database2[!database2$Group.1 %in% varRemove$chrPos,]
variantInfo = data.frame(SNP = database3$Group.1,  GENE=database3$GENE, MAF = database3$MAF_EXAC)

## raw file
datad <- read.table(paste(extract_output, ".raw", sep=""), header=T, colClasses=c(rep("character", 4), "integer", "numeric", rep("integer", nrow(mapInfo))))
### replace the phenotype with the phenotype in the pheno_file
phenoInfo = read.table(pheno_file, header=T)
phenoNew = match(phenoInfo$IID, datad$IID)
datad$PHENOTYPE[phenoNew] = phenoInfo[,eval(pheno)]-1

#### third compute the GESE test
if(familyWeightF!="None")
{	## we need to read in the weighting scheme
	familyWeight = read.table(familyWeightF, header=T, check.name=FALSE)
	gese <- GESE(pednew,  variantInfo, datad, mapInfo, threshold=1e-7, onlySeg=F, familyWeight=familyWeight )
	geseOutput <- paste(strsplit(familyWeightF, "[.]")[[1]][1], "_gese_weighted.csv", sep="")
	write.csv(gese$results[order(gese$results$pvalue_weighted),], geseOutput)

}else
{

	if(segOnly)
	{	if(recessive)
		{
			segInfo = getSegInfo(pednew, datad, mapInfo)
			### only write the genes and variants with more than 0 segregating families
			write.csv(segInfo$geneSeg[segInfo$geneSeg$numSegFam>0,], paste(extract_output, "_recessive_geneSeg.csv", sep=""))
			write.csv(segInfo$varSeg[segInfo$varSeg$numSegFam>0,], paste(extract_output, "_recessive_varSeg.csv", sep=""))
		}
	
		if(CH)
		{
			segInfo2 = getSegInfo(pednew, datad, mapInfo, mode="CH")
			### only write the genes and variants with more than 0 segregating families
			if(!anyNA(segInfo2$geneSeg))
			{	write.csv(segInfo2$geneSeg[segInfo2$geneSeg$numSegFam>0,], paste(extract_output, "_CH_geneSeg.csv", sep=""))
			}
			if(!anyNA(segInfo2$genePairSeg))
			{	write.csv(segInfo2$genePairSeg[segInfo2$genePairSeg$numSegFam>0,], paste(extract_output, "_CH_genePairSeg.csv", sep=""))
			}
			write.csv(segInfo2$varSeg[segInfo2$varSeg$numSegFam>0,], paste(extract_output, "_CH_varSeg.csv", sep=""))
		}
		
		if(dominant)
		{		
			segInfo = getSegInfo(pednew, datad, mapInfo, mode="dominant")
			### only write the genes and variants with more than 0 segregating families
			write.csv(segInfo$geneSeg[segInfo$geneSeg$numSegFam>0,], paste(extract_output, "_dominant_geneSeg.csv", sep=""))
			write.csv(segInfo$varSeg[segInfo$varSeg$numSegFam>0,], paste(extract_output, "_dominant_varSeg.csv", sep=""))

		}
	}else
	{
		gese <- GESE(pednew,  variantInfo, datad, mapInfo, threshold=1e-7, onlySeg=F)
		write.csv(gese$segregation, paste(extract_output, "_geneSeg.csv", sep=""))
		write.csv(gese$varSeg, paste(extract_output, "_varSeg.csv", sep=""))
		write.csv(gese$results[order(gese$results$pvalue),], paste(extract_output, "_gese.csv", sep=""))
	}
	
}

