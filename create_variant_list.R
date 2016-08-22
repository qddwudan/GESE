### This script creates a plink range file and a recode allele file
### input
## @ an info file for the filtered variants in the study
## @ file name for the plink range file
## @ file name for the recode allele file

args<-commandArgs(TRUE)
options(stringsAsFactors=FALSE)

outputfile = args[1]
idfile = args[2]
bimfile = args[3]
allelefile = args[4] 
snpInfofile= args[5]
freq_output2= args[6]
control_removeF= args[7]
maf_max_control= args[8]
databaseAnno = args[9]

varFile <- read.table(outputfile, sep="\t", header=T)
varFile$CHR2 = ifelse(varFile$CHR=="X", 23, ifelse(varFile$CHR=="Y", 24, varFile$CHR))

## recode-allele file
### SNP  allele
bim <- read.table(bimfile)
bim$chrPos <- paste(bim$V1, bim$V4, sep=":")
cat("****The total number of variants in the original file is ", nrow(bim), "\n")
varFile$chrPos <- paste(varFile$CHR2, varFile$POS, sep=":")
bim2 <- merge(bim, varFile, by="chrPos", sort=F)
bim3 = bim2[(bim2$ALT==bim2$V5|bim2$ALT==bim2$V6)&bim2$GENE!=".",]
cat("****The number of variants after filtering on annotation information is ", nrow(bim3), "\n")



### frequencity output file
if(file.exists(freq_output2)) 
{    ## if we were able to obtain the frequency information, possibly no control in the data to obtain this information
	freq = read.table(freq_output2, header=T)
	bim4 <- merge(bim3, freq, by.x="V2", by.y="SNP", sort=F, all.x=T)
	bim4$ALT_MAF <- ifelse(bim4$A1==bim4$ALT, bim4$MAF, ifelse(bim4$A2==bim4$ALT, 1-bim4$MAF, NA))
	maf_control = maf_max_control
	bim5 <- bim4[bim4$ALT_MAF < as.numeric(maf_control),] ## use 1% as cutoff here by default
	cat("****Remove ", sum(bim4$ALT_MAF >= as.numeric(maf_control)), " variants based on study control MAF cutoff ", maf_control, "\n")
	cat("****The number of variants after removing variants with high frequency in control is ", nrow(bim5), "\n")
	### variants removed using study control
	bim55 = bim4[bim4$ALT_MAF >= as.numeric(maf_control),]
	write.table(bim55[,c("V2", "chrPos", "CHR2", "POS", "REF", "ALT", "MAF_EXAC", "GENE", "ALT_MAF", "NCHROBS")], control_removeF, quote=F, row.names=F, col.names=T)
	

}else
{	bim5 = bim3
	bim55 = c("V2", "chrPos", "CHR2", "POS", "REF", "ALT", "MAF_EXAC", "GENE", "ALT_MAF", "NCHROBS")
	write.table(t(bim55), control_removeF, quote=F, row.names=F, col.names=T)
}

##### we want only the autosomes
bim6 <- bim5[!bim5$V1 %in% c(23, 24),]
cat("****The number of variants after removing chr 23, 24 is ", nrow(bim6), "\n")
database <- read.table(databaseAnno, sep="\t", header=T)
### If no variant present for the gene from our study in the database, likely to be coverage problem, remove these genes from our study
allGene = unique(bim6$GENE)
geneMissing = allGene[!allGene %in% database$GENE]
bim7 = bim6[!bim6$GENE %in% geneMissing,]

### range file format
## CHR START_POS END_POS SETNAME
#rangeF <- bim7[,c("CHR2", "POS", "POS", "GENE")]
rangeF <- bim7[,c("V2")]
write.table(rangeF, idfile, quote=F, row.names=F, col.names=F)

allele <- bim7[,c("V2", "ALT")]
if(file.exists(freq_output2))
{  snpInfo <- bim7[,c("V2", "GENE", "chrPos", "CHR2", "POS", "REF", "ALT", "MAF_EXAC", "ALT_MAF", "NCHROBS")]
}else
{	snpInfo <- bim7[,c("V2", "GENE", "chrPos", "CHR2", "POS", "REF", "ALT", "MAF_EXAC")]
}
colnames(snpInfo)[1:2] <- c("SNP", "GENE")
snpInfo <- snpInfo[!duplicated(snpInfo$SNP),]
write.table(allele, allelefile, quote=F, row.names=F, col.names=F)
write.table(snpInfo, snpInfofile, quote=F, row.names=F, col.names=T)

