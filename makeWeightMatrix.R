
################################################
#### A function to create the weight matrix
##############################################
createWeightMatrix <- function(extract_output=NA, studyAnno=NA, databaseAnno=NA, geneSeg_file=NA, varSeg_file=NA, pheno_file, pheno, qpheno)
{	if(anyNA(pheno_file)|anyNA(qpheno)|anyNA(pheno))
	{	stop("The phenotype file and the name of the column for the quantitative phenotype is required!")
	}
	geneSeg = read.csv(geneSeg_file, check.names=FALSE)
	fams <- colnames(geneSeg)[-c(1,2, ncol(geneSeg))]
	nfam = length(fams)
	subjInfo = read.table(pheno_file, header=T)
		
	if((anno_name)=="None")
	{	
		## the weights are the same for all the genes
		familyWeight = rep(0, nfam)
		for(i in 1:nfam)
		{	cases = subjInfo[subjInfo$FID==fams[i]&subjInfo[,eval(pheno)]==2&!is.na(subjInfo[,eval(pheno)]),]
			familyWeight[i] = mean(cases[,eval(qpheno)], na.rm=TRUE)
		}
		familyWeight2 = familyWeight/sum(familyWeight)
		familyWeight3 = data.frame(FID=fams, weight=familyWeight2)
		return (list(weight=familyWeight3))
	
	}
	
}

args <- commandArgs(TRUE)
extract_output = args[1]
studyAnno = args[2]
databaseAnno= args[3]
geneSeg_file=args[4]
varSeg_file= args[5]
pheno_file = args[6]
pheno = args[7]
qpheno = args[8]
weightF = args[9]
options(stringsAsFactors=FALSE)
resultM = createWeightMatrix(extract_output, studyAnno, databaseAnno, geneSeg_file, varSeg_file, pheno_file, pheno, qpheno)
write.table(resultM$weight, weightF, quote=F, row.names=F, col.names=T)

