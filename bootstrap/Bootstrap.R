##################################################################
## Chris Harvey
## Version 3
## February 13, 2013
##################################################################
setwd("~/Documents/classes.2012.2013/winter2013/bst.699/699.project.2/proj.2.code")
source("project.2.functions.R")
##################################################################
# global variables
g.l<-c("351", "966", "2979", "6224", "8462", "9014", "9927", "10232", "27341", "51122", "55253", "57050", "79568", "80321", "84186", "84321", "91695", "121355", "152992", "441250")	# gene labels
n.genes<-length(g.l)	# total number of genes
##################################################################
# merge data sets into a list
# data structure: dex, EtOH, & SNP.dosage
for(i in 1:n.genes){
	if(i==1){all.data=list()}
	assign("temp.data", read.table(paste("~/Documents/classes.2012.2013/winter2013/bst.699/699.project.2/proj.2.data/eqtl_data/", g.l[i], ".eqtl.dat", sep=""), header=TRUE))
	all.data[[i]]<-temp.data
	names(all.data)[[i]]<-paste("data.", g.l[i], sep="")
}
##################################################################
## reliant on the list of all data 'all.data'
for(i in 1:n.genes){
	# i<-1
	
	if(i==1){all.tTables=list()}		# list of tTables for every gene
	
	n.subjects<-dim(all.data[[i]])[1]		# number of subjects per gene (should always be 57)
	n.snps<-dim(all.data[[i]])[2]-2			# number of SNPs in current gene
	B.bootstraps<-1050						# number of bootstrap samples
	
	# reshape data into long format and add an indicator for treatment
	data.orig.long<-reshape.long.treat(all.data[[i]], n.subjects)
	
	# data structure for tTable with original model estimates and bootstrap SE, t-val, and P
	orig.boot.tTable<-as.data.frame(matrix(nrow=4*n.snps, ncol=7))
	orig.boot.tTable<-cbind(rep(c("B.o", "B.SNP", "B.dex", "B.SNP.dex"), 4*n.snps), orig.boot.tTable)
	names(orig.boot.tTable)<-c("predictor", "estimate", "SE", "t-val", "P(>|t|)", "bs.SE", "bs.t-val", "bs.P(>|t|)")
	
	print(Sys.time())
	for(j in 1:n.snps){
		# j<-1
		print(paste("SNP = ", j))
		
		# increment row multiplier in the next line
		if(j==1){index<-0}		
		rows<-c((index*4+1):(index*4+4))
		
		#print(rows)
		
		# full model from original sample(expression ~ currentSNP + dex + currentSNP*dex)
		original.SNP.model<-lm(data.orig.long[, 1]~ data.orig.long[, j+2] + data.orig.long[, 2] + data.orig.long[, j+2]*data.orig.long[, 2])
		orig.boot.tTable[rows,2:5]<-summary(original.SNP.model)$coefficients
		
		for(b in 1:B.bootstraps){
			# b<-3
			
			 # initialize matrice for bootstrapped coefficients (intercept, B.SNP, B.dex, B.SNP.dex)      
			if(b==1){ boot.coefficients<-as.data.frame(matrix(nrow = B.bootstraps, ncol = 4))}	
			boot.sample<-as.integer(sample(seq(1:n.subjects), n.subjects, replace=TRUE))		# bootstrap of original sample
			temp.boot.data<-all.data[[i]][boot.sample, c(1,2,j+2)]								# derive data from bootstrap and original data
			temp.boot.data.long<-reshape.long.treat(temp.boot.data, n.subjects)					# shape new data into long format, add a treatment term
			
			# model from boot-strapped sample(expression ~ currentSNP + dex + currentSNP*dex)
			boot.model<-lm(temp.boot.data.long[, 1]~ temp.boot.data.long[, 3] + temp.boot.data.long[, 2] + temp.boot.data.long[, 3]*temp.boot.data.long[, 2])
			
			# store coefficients
			boot.coefficients[b, ]<-boot.model$coefficients
			
		} # B.bootstraps
		
		orig.boot.tTable[rows,6]<-sapply(na.omit(boot.coefficients), sd)
		orig.boot.tTable[rows,7]<-orig.boot.tTable[rows,2]/orig.boot.tTable[rows,6]
		orig.boot.tTable[rows,8]<-t.pval(orig.boot.tTable[rows,7], 110)
		
		print(orig.boot.tTable[rows,1:8])
		
		index<-index+1
	
	} # n.snps
	
	print(Sys.time())
	all.tTables[[i]]<-orig.boot.tTable
	write.csv(orig.boot.tTable, paste("BS.tTable.gene.", g.l[i], ".csv", sep=""))
	
	if(i==n.genes){names(all.tTables)=g.l}
	
} # n.genes

