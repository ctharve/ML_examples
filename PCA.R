##################################################################
## Chris Harvey
## Version 2
## April 13, 2013
##################################################################
setwd("~/Documents/classes.2012.2013/winter2013/bst.699/699.project.5/proj.5.code")
#################################
data.sample<-read.table("~/Documents/classes.2012.2013/winter2013/bst.699/699.project.5/proj.5.data/sample_info.dat", header=TRUE)
data.otu<-read.table("~/Documents/classes.2012.2013/winter2013/bst.699/699.project.5/proj.5.data/OTUs.shared.dat", header=TRUE)
data.mice<-read.table("~/Documents/classes.2012.2013/winter2013/bst.699/699.project.5/proj.5.data/mice_weights.dat", header=TRUE)
##################################################################
KrusPVal<-function(pc, indic){
	return(kruskal.test(pc, indic)$p.val)	
}
##################################################################
library(lme4)
library(gplots)
##library(coefplot2)		## this package is not available through CRAN
library(ggplot2)
library(reshape)
library(RColorBrewer)
library(scatterplot3d)
#################################
setwd("~/Documents/classes.2012.2013/winter2013/bst.699/699.project.4/project.4.output")
##################################################################
head(data.otu)
data.otu[1:3, 1:3]
dim(data.otu)
matrix.design<-data.otu[, -c(1,2)]
n<-dim(data.otu)[1]
p<-dim(matrix.design)[2]
s.hat<-cov(data.otu[, -c(1,2)])
s.hat.2<-s.hat%*%s.hat
s.hat.trace<-sum(diag(s.hat))
s.hat.2.trace<-sum(diag(s.hat.2))
v.hat<-s.hat.trace/p
f.hat<-v.hat*diag(p)

row.test<-((1-2/p)*s.hat.2.trace+s.hat.trace^2)/((n+1-2/p)*(s.hat.2.trace-s.hat.trace^2/p))
row.hat<-min(row.test, 1)

s.hat.oas<-row.hat*f.hat+(1-row.hat)*s.hat
print(s.hat.oas)
dim(s.hat.oas)

s.oas.eigen<-eigen(s.hat.oas)
dim(s.oas.eigen$vectors)

## only 15 eigenvalues because n=15
s.oas.eigen$values

eigen.values<-s.oas.eigen$values[1:15]
eigen.porportion<-eigen.values/sum(eigen.values)

cummulative.eigen.sums<-cumsum(eigen.porportion)

pc.1<-as.matrix(matrix.design)%*%s.oas.eigen$vectors[, 1]
pc.2<-as.matrix(matrix.design)%*%s.oas.eigen$vectors[, 2]
pc.3<-as.matrix(matrix.design)%*%s.oas.eigen$vectors[, 3]
pc.4<-as.matrix(matrix.design)%*%s.oas.eigen$vectors[, 4]
pc.5<-as.matrix(matrix.design)%*%s.oas.eigen$vectors[, 5]
pc.6<-as.matrix(matrix.design)%*%s.oas.eigen$vectors[, 6]

### make 3-D R plot of the first three pcs and label the points, size of points differs based on exposure conditions
pcs<-as.data.frame(cbind(pc.1, pc.2, pc.3, pc.4, pc.5, pc.6))
names(pcs)<-c("ev1", "ev2", "ev3", "ev4", "ev5", "ev6")

attach(pcs)
cube<-scatterplot3d(pcs$ev1, pcs$ev2, pcs$pcsev3)

detach(pcs)

##################################################################
data.otu[, 1:2]
data.sample
## 381lg missing OTU data
data.sample.sort<-data.sample[order(data.sample[, 1]), ]
exp.ppm.ind<-c(1, 1, 1, 3, 3, 3, 3, 1, 1, 1, 3, 2, 2, 3, 1)
length(exp.ppm.ind)

## testing for seggregation based on each OTU
p.values.nonPara<-apply(matrix.design, 2, KrusPVal, indic=exp.ppm.ind)
hist(p.values.nonPara)
range(p.values.nonPara)
p.vales.bonferonni<-(p.values.nonPara<.05/length(p.values.nonPara))

##################################################################
head(data.mice)
qplot(week, weight, data=data.mice, colour=id)

data

data.mice<-cbind(data.mice, rep(NA, 228))
names(data.mice)[4]<-"weight.adjust"
data.mice[1:34, 4]<-data.mice[1:34, 3]/7.35
data.mice[35:64, 4]<-data.mice[35:64, 3]/7.19
data.mice[65:93, 4]<-data.mice[65:93, 3]/7.68
data.mice[94:123, 4]<-data.mice[94:123, 3]/7.47
data.mice[124:158, 4]<-data.mice[124:158, 3]/9.2
data.mice[159:192, 4]<-data.mice[159:192, 3]/8.54
data.mice[193:228, 4]<-data.mice[193:228, 3]/9.40

qplot(week, weight.adjust, data=data.mice, colour=id)
attach(data.mice)
data.mice<-data.mice[order(id, week), ]
detach(data.mice)

levels(data.mice$id)

data.mice<-cbind(seq(1, dim(data.mice)[1]), data.mice)

dim(data.mice)
data.mice<-cbind(data.mice, c(rep(0,94), rep(1,134)))
names(data.mice)[6]<-"exposure"
qplot(week, weight.adjust, data=data.mice, colour=exposure)
data.mice$exposure<-as.factor(data.mice$exposure)
p.tmp<-ggplot(data.mice, aes(week, weight.adjust, colour=id, lty=exposure))
p.tmp+layer(geom="line")

##################################################################
