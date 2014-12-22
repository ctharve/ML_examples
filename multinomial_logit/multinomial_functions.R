###########################################################################
## multinomial_functions.R
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
## Functions:
## + YHatMulti (u_ic in Murphy's MLAPP)
## + MultiNegLlk
## + GradMultiNegLlk
## + HessMultiNegLlk
##
###########################################################################

## YHatMulti  
## w ~ @param length P*(C-1) vector of parameters 
## X ~ @param N X P matrix of obervations
## baseline ~ @param boolean, include prediction of baseline class
## Y_hat ~ @val N X C matrix of predicted class probabilities
YHatMulti <- function(w, X, baseline=TRUE){
  p <- dim(X)[2]
  W <- matrix(w, nrow=p) 
  class_inner_prod <- exp(X %*% W)
  denom <- (1 + rowSums(class_inner_prod))^-1
  if(baseline){
    cbind(class_inner_prod * denom, denom)
  } else {
    class_inner_prod * denom
  }
}

## MultiNegLlk
## w ~ @param length P*(C-1) vector of parameters 
## X ~ @param N X P matrix of obervations
## Y ~ @param N X C matrix of bits indicating class for each observation
## -llk ~ @val negative log-likelihood evaluated at w 
MultiNegLlk <- function(w, X, Y){
  Y_hat <- YHatMulti(w, X)  
  llk <- sum(diag(t(Y) %*% log(Y_hat)))  # trace is invariant to cyclic permutations
  -llk
}

## GradMultiNegLlk
## w ~ @param length P*(C-1) vector of parameters 
## X ~ @param N X P matrix of obervations
## Y ~ @param N X C matrix of bits indicating class for each observation
## grad ~ @val gradient vector of length P*(C-1)   
GradMultiNegLlk <- function(w, X, Y){
  n <- dim(Y)[1]
  c <- dim(Y)[2]
  p <- dim(X)[2]
  if(n != dim(X)[1]){stop("Err: inconsistent input")}
  Y_hat <- YHatMulti(w, X)
  Y_error <- (Y_hat - Y)[, -c] ## (N X C-1) ## category C is baseline
  k_prod <- sapply(1:n, function(ii){
    kronecker(Y_error[ii, ], X[ii, ])
  })
  if(dim(k_prod)[1] != (p * (c - 1))){stop("Kronecker product has inconsitent row dimension")}
  if(dim(k_prod)[2] != n){stop("Kronecker product has inconsitent column dimension")}
  grad <- rowSums(k_prod)
  -grad
}

## WARNING: this is buggy
## Multinomial Model: hessian of the log-likelihood
## w ~ @param length P*C vector of parameters (include class specific intercepts)
## X ~ @param N X P matrix of obervations
## hess ~ @val P*(C-1) hessian    
HessMultiNegLlk <- function(w, X){
  n <- dim(X)[1]
  p <- dim(X)[2]
  Y_hat <- YHatMulti(w, X, baseline=FALSE)  ## N X C-1 matrix of predictions
  hess_list <- vector("list", n)

  ## this might take forever. like, actually forever. 
  invisible(sapply(1:n, function(ii){
    Y_hat_vec  <- Y_hat[ii, ]
    X_vec <- X[ii, ]
    X_outer_prod <- X_vec %*% t(X_vec)
    ##DANGER## k_prod <- kronecker(diag(Y_hat_vec), X_outer_prod)
    hess_list[[ii]] <<- k_prod
  }))
  
  hess <- Reduce("+", hess_list)
  -hess
}

##         ##
## THE END ##
##         ##
