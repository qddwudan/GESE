#### This script is used to create subject FID IID file for plink to extract

### input
## pheno_file name
## pheno: phenotype column name
## output_file name

args <- commandArgs(TRUE)
options(stringsAsFactors=FALSE)
pheno_file = args[1]
pheno = args[2]
output_file=args[3]
output_file_control=args[4]

phenoInfo <- read.table(pheno_file, header=T)
pheno2 = phenoInfo[phenoInfo[,eval(pheno)] %in% c(1,2), 1:2]
write.table(pheno2, output_file, quote=F, row.names=F, col.names=F)

### controls 
controls = phenoInfo[phenoInfo[, eval(pheno)] == 1, 1:2]
if(nrow(controls)>0)
{
write.table(controls, output_file_control, quote=F, row.names=F, col.names=F)
}
