###########################################################################
## multinomial_testing.R
## Multinomial logit model ~ negative log-likelihood, gradient, and hessian
##
## CTH 12.22.14
##
## Reference:
## Machine Learning: A Probabilistic Perspecitve by Kevin P. Murphy 2012
##
## About:
## logit (multinomial) regression model with
## N ~ individuals
## C ~ categories
## P ~ parameters (including intercept)
##
## Analysis:
## testing multinomial_functions.R for numerical stability
## testing multinomial_optimization.R for run-time. Potential to parallelize by
## partitioning the data and running for a few iterations, then restarting all
## partitions on estimates that maximize the total likelihood. RPR 121914
##
###########################################################################
library(VGAM)
source("./MESH_multinomial_functions.R")
#X11(type="Xlib", display="localhost:12.0")

## Simulate the data 
## http://www.uni-kiel.de/psychologie/rexrepos/posts/regressionMultinom.html

## test data
set.seed(666)
N  <- 100
X1 <- rnorm(N, 175, 7)
X2  <- rnorm(N,  30, 8)
Ycont <- 0.5*X1 - 0.3*X2 + 10 + rnorm(N, 0, 4)
Ycateg <- cut(Ycont, breaks=quantile(Ycont), include.lowest=TRUE, labels=c("--", "-", "+", "++"), ordered=FALSE)
dfMN <- data.frame(X1, X2, Ycateg)

## model fit using vglm
vglmFitMN <- vglm(Ycateg ~ X1 + X2, family=multinomial(refLevel=1), data=dfMN)  
exp(VGAM::coef(vglmFitMN))


## test data for optimization with optim
testDat <- as.matrix(dfMN[, -3]) ## remove the labels 
testDat <- cbind(rep(1, dim(testDat)[1]), testDat)
outcome <- as.character(dfMN$Ycateg)
Y <- Reduce(rbind, lapply(1:length(outcome), function(ii){
   switch(outcome[ii],
          "--" = c(0, 0, 0, 1),
          "-" = c(1, 0, 0, 0),
          "+" = c(0, 1, 0, 0),
          "++" = c(0, 0, 1, 0))
}))

p <- dim(testDat)[2]
c <- dim(Y)[2] - 1
w <- rep(0, p*c)
a
## 0.) do the functions return sane (relatively) results.
myllk <- MultiNegLlk(w, testDat, Y)
mygrad <- GradMultiNegLlk(w, testDat, Y)  
myhess <- HessMultiNegLlk(w, testDat)
myllk; mygrad; myhess

## 1.) optimization: llk
testoptim_llk <- optim(w, fn=MultiNegLlk, X=testDat, Y=Y,method="BFGS")
otherpar <- VGAM::coef(vglmFitMN)[c(1,4,7,2,5,8,3,6,9)]
testoptim_llk$par; otherpar
cat('my nllk:', MultiNegLlk(testoptim_llk$par, testDat, Y), '\n')  ## 125.3079
cat('vgam nllk:', MultiNegLlk(otherpar, testDat, Y), '\n')  ## 124.9729
cat('nllk at 0', MultiNegLlk(w, testDat, Y), '\n')  ## 138.629
cat('nllk at 1', MultiNegLlk(rep(1, p*c), testDat, Y), '\n')  ## 5253.439

## 2.) optimization: llk + g_llk
testoptim <- optim(w, fn=MultiNegLlk, gr=GradMultiNegLlk, X=testDat, Y=Y, method="L-BFGS") ## may be buggy
cat('my grad-nllk:', MultiNegLlk(testoptim$par, testDat, Y), '\n')  ## nllk 125.3079  
cat('vgam nllk:', MultiNegLlk(otherpar, testDat, Y), '\n')  ##nllk 124.9729

## 2.a) performance on a larger data set
sample.ind <- sample(1:n, 500000, replace=TRUE)
testDat.big <- testDat[sample.ind, ]
Y.big <- Y[sample.ind, ]
s.time <- proc.time()
testoptim.big <- optim(w, fn=MultiNegLlk, gr=GradMultiNegLlk, X=testDat.big, Y=Y.big, method="L-BFGS") ## may be buggy
e.time <- proc.time()
cat(e.time-s.time); # 802.355 seconds
cat('my grad-nllk:', MultiNegLlk(testoptim.big$par, testDat.big, Y.big), '\n')  ## nllk 626443.3

dfMN.big <- dfMN[sample.ind, ]
s.time <- proc.time()
vglmFitMN.big <- vglm(Ycateg ~ X1 + X2, family=multinomial(refLevel=1), data=dfMN.big)  
e.time <- proc.time()
cat(e.time-s.time); #32 seconds
otherpar.big <- VGAM::coef(vglmFitMN.big)[c(1,4,7,2,5,8,3,6,9)]
cat('vgam nllk:', MultiNegLlk(otherpar.big, testDat.big, Y.big), '\n')  ##nllk 624785.4


## 3.) gradient diagnostics
wdummy<- rep(0, 9)
wtest <- testoptim_llk$par ## parameters from the llk only model
nllk <- MultiNegLlk(wtest, testDat, Y)
gtest <- GradMultiNegLlk(wtest, testDat, Y)
gamma <- .1
index <- 1

w <- wtest
X <- testDat




## plot the vector field
index = 4
gamma <- (-100:100)/40
nllk_pert <- sapply(gamma,function(gamma) MultiNegLlk(MutVec(wtest, index, gamma), testDat, Y))
gr_pert <- sapply(gamma,function(gamma) GradMultiNegLlk(MutVec(wtest, index, gamma), testDat, Y))
x_pert <- sapply(gamma,function(gamma) MutVec(wtest, index, gamma))
x_inc <- sapply(gamma,function(gamma) MutVec(wtest, index, gamma) + GradMultiNegLlk(MutVec(wtest, index, gamma), testDat, Y)*0.00001)
y_inc <- sapply(gamma,function(gamma) sum(GradMultiNegLlk(MutVec(wtest, index, gamma), testDat, Y)^2)*0.0001^2)
ii <- 1:length(gamma)
plot(x_pert[index,],nllk_pert,type='l')
arrows(x_pert[index,ii],nllk_pert[ii], x_inc[index,ii],nllk_pert[ii],length=0.1)

plot(gamma,gr_pert[index,],type='l')

cat('original llk: ', nllk, '\n');
cat('perturbed llk: ', nllk_pert, '\n');
cat('finite diff llk: ', nllk + gtest%*%MutVec(wdummy, index, gamma), '\n')

MutVec <- function(x, index, gamma){
  x[index] = x[index] + gamma
  x
}



## 3.) compare model fit predictions: use parameters to predict outcomes on test data 
testoptim_llk_preds <- YHatMulti(testoptim_llk$par, testDat)
testoptim_preds <- YHatMulti(testoptim$par, testDat)
other_preds <- YHatMulti(otherpar, testDat)

##
ps2 <- VGAM::predict(vglmFitMN)
mu1 <- 1/(1+rowSums(exp(ps2)))
preds2 = cbind(exp(ps2)*mu1,mu1)
t(preds2) %*% Y



##         ##
## THE END ##
##         ##
