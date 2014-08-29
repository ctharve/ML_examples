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
hist(x[, "cylinders"], breaks=10, main='cylinders')  ## discrete [3, 8]
hist(x[, "displacement"], breaks=100, main='displacement') ## cont [~45, 400+]
hist(x[, "horsepower"], breaks=100, main='horsepower') ## cont [~48, 200+]
hist(x[, "weight"], breaks=50, main='weight')  ## cont [1500, 5000]
hist(x[, "acceleration"], breaks=50, main='acceleration')  ## cont [~0, 25]
hist(y, breaks=50, main='mpg')  ## cont [5, 45]

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
cost <- vector(mode='integer', B)   ## monitor the cost/objective function
colnames(thetas) <- colnames(xn)   ## theta for each column in data matrix
final_iter <- NA

##
for(tt in 2:B){
  #tt <- 2
  prev <- tt-1
  #cat(tt, prev, '\n')
  
  ## will clean this up after testing
  y_hat <- xn %*% thetas[prev, ]
  diff <- (y_hat - y)
  cost[tt] <- n.obs^(-1)*sum(diff^2)
  thetas[tt, ] <- thetas[prev, ] - alpha*(n.obs^(-1))*( t(xn) %*% diff)
  #thetas[tt, 'cylinders'] <- thetas[prev, 'cylinders'] - alpha*(n.obs^(-1))*( t(x[, 'cylinders']) %*% (y_hat - y))    
  
  ## reporting
  if(tt%%100 == 0){cat(paste('iteration', tt, sep=' '), '\n', colnames(thetas), '\n', thetas[tt, ], '\n', 'Cost function=', cost[tt], '\n')}
  
  ## convergence
  if(sum( abs(thetas[tt, ]-thetas[prev, ]) <= epsilon) == n.thetas){
    final_iter <- tt
    cat('Converged at iteration ', tt)
    break
  }

}

thetas[c((final_iter-5):final_iter), ]
plot(cost[c(2:100)], main='cost')
plot(thetas[c(2:400), 'intercept'], main='intercept')
plot(thetas[c(2:400), 'cylinders'], main='cylinders')
plot(thetas[c(2:400), 'displacement'], main='displacement')
plot(thetas[c(2:400), 'horsepower'], main='horsepower')
plot(thetas[c(2:400), 'weight'], main='weight')
plot(thetas[c(2:400), 'acceleration'], main='acceleration')

##               ##
## ## THE END ## ##
##               ##  