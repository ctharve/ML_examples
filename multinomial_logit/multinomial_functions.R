###########################################################################
##
## CTH 120814
## Multinomial logit model
##
## functions for a generalized logit (multinomial) model with
## N ~ individuals
## C ~ categories
## P ~ parameters including intercept
##
## + Y_hat (u_ic in Murphy's MLAPP)
## + negative log-likelood
## + gradient
## + hessian
##
###########################################################################
## YHatMulti  
## w ~ @param length P*(C-1) vector of parameters 
## X ~ @param N X P matrix of obervations
## baseline ~ @param boolean, include prediction of baseline class
## Y_hat ~ @val N X C matrix of predicted class probabilities
YHatMulti <- function(X, w, baseline=TRUE){
  p <- dim(X)[2]
  W <- matrix(w, nrow=p)
  
  class_inner_prod <- exp(X %*% W)
  denom <- (1 + rowSums(class_inner_prod))^-1
  
  if(baseline){
    Y_hat <- cbind(class_inner_prod * denom, denom)
    Y_hat
  } else {
    Y_hat <- class_inner_prod * denom
    Y_hat
  }
}

## MultiNegLlk
## w ~ @param length P*(C-1) vector of parameters 
## Y ~ @param N X C matrix of bits indicating class for each observation
## X ~ @param N X P matrix of obervations
## -llk ~ @val negative log-likelihood evaluated at w 
MultiNegLlk <- function(w, X, Y){
  Y_hat <- YHatMulti(X, w)  
  llk <- sum(log(diag(t(Y) %*% Y_hat)))  # trace is invariant to cyclic permutations
  -llk
}

## GradMultiNegLlk
## w ~ @param length P*(C-1) vector of parameters 
## Y ~ @param N X C matrix of bits indicating class for each observation
## X ~ @param N X P matrix of obervations
## grad ~ @val gradient vector of length P*(C-1)   
GradMultiNegLlk <- function(w, X, Y){
  n <- dim(Y)[1]
  c <- dim(Y)[2]
  p <- dim(X)[2]
  if(n != dim(X)[1]){stop("Err: inconsistent input")}
  Y_hat <- YHatMulti(X, w)
  Y_error <- t((Y_hat - Y)[, -c]) ## (C-1 X N) ## C^{th} category is baseline
  k_prod <- kronecker(Y_error, t(X))
  if(dim(k_prod)[1] != (p * (c - 1)){stop("Err: kronecker product has inconsitent dimension")}
  grad <- rowSums(k_prod)
  grad
}

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
    k_prod <- kronecker(diag(Y_hat_vec), X_outer_prod)
    hess_list[[ii]] <<- k_prod
  }))
  
  hess <- Reduce("+", hess_list)
  hess
}

##         ##
## THE END ##
##         ##
