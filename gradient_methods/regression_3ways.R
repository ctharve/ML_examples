##########################################################################################
## protoyping gradient based methods before implementation in Java
##########################################################################################
library(dplyr)
library(ggplot2)

## load data 
##  
dat <- read.table('Desktop/sgd_fun/sgd_R/auto-mpg.data.txt', stringsAsFactors=FALSE)
names(dat) <- c('mpg', 'cylinders', 'displacement', 'horsepower', 'weight', 'acceleration', 'model-year', 'origin', 'car-name')
dat_clean <- dat[-which(dat$horsepower=='?'), ]
dat_clean$horsepower <- as.numeric(dat_clean$horsepower)
dat_clean$cylinders <- as.numeric(dat_clean$cylinders)
head(dat_clean)

## 1.) learning with lm package
##
pack_model <- lm(mpg ~ cylinders + displacement + horsepower + weight + acceleration, dat_clean)
summary(pack_model)
coefs_pack <- pack_model$coefficients; coefs_pack

## 2.) learning with normal equations
##
drop <- c('mpg', 'model-year', 'origin', 'car-name')
feature_names <- names(dat)[!(names(dat_clean)%in% drop)]
intercept <- as.numeric(rep(1, dim(dat_clean)[1]))
x <- as.matrix(cbind(intercept, dat_clean[, feature_names])) 
head(x)
str(x)
## for gradient descent, it's a good idea to mean center and scale all of the parameters, 
## but first, let's do what we should have done in the first place, which is to plot all of the features.
##
hist(x$cylinders, breaks=10, main='cylinders')  ## discrete [3, 8]
ggplot(x, aes(x=cylinders)) + geom_density()
hist(x$horsepower, breaks=100, main='horsepower') ## cont [~48, 200+]
hist(x$displacement, breaks=100, main='displacement') ## cont [~45, 400+]
hist(x$weight, breaks=50, main='weight')  ## cont [1500, 5000]
## xn ~ normalized features
xn <-  cbind(intercept, Reduce(cbind, lapply(c(2:6), FUN=function(this_feature){
    mu <- mean(x[, this_feature])
    se <- sd(x[, this_feature])
    '*'('-'(x[, this_feature], mu), (se^-1)) 
    })
  )
)

features_full <- c('intercept', feature_names)
colnames(x) <- features_full
colnames(xn) <- features_full
y <- dat_clean[, 'mpg']
beta <- (solve(t(x) %*% x))%*%(t(x)%*%y); beta 
beta_n <- (solve(t(xn) %*% xn))%*%(t(xn)%*%y); beta_n 

## 3.) learning with gradient descent 
##
n.obs <- dim(x)[1]    ## number of observations
n.thetas <- dim(x)[2]   ## number of parameters
B <- 5000    ## max iterations
epsilon <- 1E-6   ## convergence threshold
alpha <- 0.05    ## learning rate
thetas <- matrix(0, B, n.thetas)   ## initialize matrix for results
thetas[1, ] <- rep(1, n.thetas)     ## chose zeroes as starting points (may want to change this)
#thetas[1, ] <- beta
colnames(thetas) <- colnames(x)   ## theta for each column in data matrix

# replace the original data with normmalized data
x_original <- x
x <- xn
#

for(tt in 2:B){
  #tt <- 2
  prev <- tt-1
  #cat(tt, prev, '\n')
  
  ## will clean this up after testing
  thetas[tt, 'intercept'] <- thetas[prev, 'intercept'] - alpha*(n.obs^(-1))*( t(x[, 'intercept']) %*% (x %*% (thetas[prev, ]) - y))    
  thetas[tt, 'cylinders'] <- thetas[prev, 'cylinders'] - alpha*(n.obs^(-1))*( t(x[, 'cylinders']) %*% (x %*% (thetas[prev, ]) - y))    
  thetas[tt, 'displacement'] <- thetas[prev, 'displacement'] - alpha*(n.obs^(-1))*( t(x[, 'displacement']) %*% (x %*% (thetas[prev, ]) - y))    
  thetas[tt, 'horsepower'] <- thetas[prev, 'horsepower'] - alpha*(n.obs^(-1))*( t(x[, 'horsepower']) %*% (x %*% (thetas[prev, ]) - y))    
  thetas[tt, 'weight'] <- thetas[prev, 'weight'] - alpha*(n.obs^(-1))*( t(x[, 'weight']) %*% (x %*% (thetas[prev, ]) - y))    
  thetas[tt, 'acceleration'] <- thetas[prev, 'acceleration'] - alpha*(n.obs^(-1))*( t(x[, 'acceleration']) %*% (x %*% (thetas[prev, ]) - y))    
  
  
  ## reporting
  ## every hundred iterations print current estimates
  if(tt%%100 == 0){cat(paste('iteration', tt, sep=' '), '\n', colnames(thetas), '\n', thetas[tt, ], '\n')}
  #if(tt%%10 == 0){cat(paste('iteration', tt, sep=' '), '\n', theta, '\n')}
  
  
  ## convergence
  if(sum( abs(thetas[tt, ]-thetas[prev, ]) <= epsilon) == n.thetas){
    cat('Converged at iteration ', tt)
    break
  }

}
head(thetas); tail(thetas)

##               ##
## ## THE END ## ##
##               ##  