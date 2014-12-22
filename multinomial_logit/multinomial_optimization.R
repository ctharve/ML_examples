###########################################################################
## multinomial_optimization.R
## Multinomial logit ~ optimization, output processing, & visualization
##
## CTH 12.22.14
## 
## About:
## logit (multinomial) regression model with
## N ~ individuals
## C ~ categories
## P ~ parameters (including intercept)
##
## Functions:
## + MungeResults
## + MyMultinomial
##
###########################################################################

## MungeResults
## optobj ~ @param an optim object from multinomial model
## parmnames ~ @param vector of all names
## classes @param number of outcome classes
## sig ~ @param level to determine significance
## retval ~ @val list with data frames of all data, data as matrices, and significant parameters
## TOD0: add forest plot of all significant data
MungeResults <- function(optobj, parmnames, classes, sig=0.05){
  n.c <- classes-1
  n.p <- length(parmnames)/n.c

  parms <- optobj$par
  parms.se <- diag(solve(optobj$hessian))^(.5)
  parms.z <- parms * parms.se^(-1)
  parms.pval <- 2*pnorm(-abs(parms.z))
  betas <- exp(parms)
  dat.parms <- data.frame(betas, parms, parms.se, parms.z, parms.pval)

  index <- which(dat.parms$parms.pval<(sig))
  dat.sig <- dat.parms[index, ]

  ## this should be done over a list, so we will fix it
  pmat.betas <- format(matrix(betas, ncol=n.c), digits=3, scientific=F)
  pmat.parms <- format(matrix(parms, ncol=n.c), digits=3, scientific=F)
  pmat.se <- format(matrix(parms.se, ncol=n.c), digits=3, scientific=F)
  pmat.z <- format(matrix(parms.z, ncol=n.c), digits=3, scientific=F)
  rownames(pmat.parms) <- parmnames[1:n.p]
  rownames(pmat.se) <- parmnames[1:n.p]
  rownames(pmat.z) <- parmnames[1:n.p]
  dat.mat <- list(pmat.betas=pmat.betas, pmat.parms=pmat.parms, pmat.z=pmat.z, pmat.se=pmat.se)
  
  retval <- list(dat.all=dat.parms, dat.mat=dat.mat, dat.sig=dat.sig)
  retval   
}

## MyMultinomial
## multnomial regression optimized on the negative-llk using BFGS
## parms ~ @param initial parameters
## X ~ @param N X P design matrix with a column of 1s for intercept
## Y ~ @param N X C matrix of outcomes (bit vectors or probabilities) 
## var.names ~ @param N*(C-1) vector of variable names (WARNING: ensure the names align with parms)
## retval ~ @val important things
##
## TODO: test with the gradient and hessian
MyMultinomial <- function(parms, X, Y, var.names, hessian=TRUE, maxit=1E5){

  ## check argument validity ##
  cat('=== parameters are valid ===', '\n\n')

  ## model fitting
  cat('=== begin optimization ===', '\n\n')
  st <- proc.time()
  model <- optim(parms, fn=MultiNegLlk, X=as.matrix(X3), Y=as.matrix(Y), method="BFGS", hessian=TRUE, control=list(maxit=maxit));
  en <- proc.time()
  cat('Run time: ', en-st, '\n\n');
  
  ## results: model fit and parameter estimates
  cat('=== processing ouput ===', '\n\n')
  nllk <- model$value 
  converged <- modle$convergence
  hessian <- model$hessian
  model_res <- MungeResults(model, var.names, classes=4)

  retval <- list(nllk=nllk, converged=converged, hessian=hessian, results=model_res)
  retval
}

##         ##
## THE END ##
##         ##
