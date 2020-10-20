# QuantileTensorReg

We present implementation for estimation of quantile tensor regression model. The estimate is obtained by an alternating update algorithm based on Tucker decomposition. We also use a sparse Tucker decomposition to further reduce the number of parameters when the dimension of the tensor is large. We propose an alternating update algorithm combined with alternating direction method
of multipliers (ADMM) for sparse scenario.

## Usage
For nonsparse case, one can use 
```{r lowfun, eval=FALSE}
lowdim_AQR(y, X, X_vec, a0, rank, tau, tol = 1e-06, iter_max)
```
Arguments: <br>
  * y:  Response vactor.
  * X: A list of input tensor.
  * X\_vec: A list of the same length as order of X[[i]]. The kth element is a $n$ row matrix of vectorized X, each row of which is obtained by vectorizing mode-k matricization of X[[i]]. X\_vec can be obtained by function X\_unfold(X,mode = k).
  * a0: Initial.
  * rank: Predetermined rank.
  * tau: The conditional quantile of the response to be estimated, so must be a number between 0 and 1.
  * tol: Convergence threshold. 
  * iter\_max: Maximum number of iterations. <br>

This function return a list containing:<br>
  * A: The estimated coefficients.
  * U: A list of estimated factor matrices.
  * G: The estimated core tensor.
  * BIC: Value of BIC.
  * iter: The number of iterations.


For sparse case, one can use 
```{r highfun, eval=FALSE}
Highdim_AQR(y, X, X_vec, a0, rank, tau, lambda, penalty = c(TRUE,TRUE,TRUE),tol = 4e-03, iter_max)
```
The two additional arguments:<br>
  * lambda: tuning parameter 
  * penalty: A vector of the same length as order of X[[i]]. Logical; if TRUE, penalty will be imposed on the corresponding factor matrix.



Some examples:

```{r example1, eval=FALSE}
###########################
### Example 1 nonsparse
###########################
library(rTensor)
library(quantreg)

tau = 0.5
dim = c(10,10,10)
rank = c(3,3,3)
n = 2000 
sigma_g = 1
sigma = 0.5
order = 3

###########generate data############
A = DGP_A(sigma_g, dim = dim, rank = rank)$A
dat = DGP(n, coef = A, sigma = sigma, tau = tau, dim = dim)
y = dat$y
X = dat$X

X_vec = list()
for(k in 1:order){
  X_vec[[k]] = X_unfold(X,mode = k)
}

##########  initial  ###############
A0 = rq.fit.lasso(X_vec[[1]],y, tau = tau, lambda = 0.0001)
a0 = array(A0$coefficients[-1],dim = dim)
a0 = as.tensor(a0)

######### select  rank  ###################
rankgrid = as.matrix( expand.grid(c(2,3,4),c(2,3,4),c(2,3,4)) )
BIC = vector()
for (i in 1:nrow(rankgrid)) {
  ranktemp = rankgrid[i,]
  try = try(lowdim_AQR(y, X, X_vec, a0, rank = ranktemp, tau = tau, tol = 1e-06, iter_max = 200))
  if("try-error" %in%class(try)){cat("error when rank =", ranktemp);next}
  BIC[i] = try$BIC
  print(i)
}

rankhat = rankgrid[which.min(BIC),]

########estimate  ################
Ahat = lowdim_AQR(y, X, X_vec, a0, rank = rankhat, tau = tau, tol = 1e-06, iter_max = 200)

fnorm(Ahat$A-A) #estimation error






###########################
### Example 2 sparse
###########################
library(rTensor)
library(quantreg)

tau = 0.5
dim = c(10,10,10)
rank = c(2,2,2)
n = 400
sigma_g = 1
sigma = 1
order = 3

tnr = DGP_Ahigh(sigma_g, dim = dim, rank = rank)
A = tnr$A
Utrue = tnr$U
Gtrue = tnr$G
dat = DGP(n, coef = A, sigma = sigma, tau = tau, dim = dim)
y = dat$y
X = dat$X

X_vec = list()
for(k in 1:order){
  X_vec[[k]] = X_unfold(X,mode = k)
}

dat_val = DGP(400, coef = A, sigma = sigma, tau = tau, dim = dim)
y_val = dat_val$y
X_val = dat_val$X
X_val_vec = X_unfold(X_val,mode = 1)

##########  initial  ###############
A0 = rq.fit.lasso(X_vec[[1]],y,lambda = 0.0001,eps = 0.1);
a0 = array(A0$coefficients[-1],dim=dim)
a0 = as.tensor(a0)

#########  select rank  ###############
rankgrid = as.matrix( expand.grid(c(1,2,3),c(1,2,3),c(1,2,3)) )

pre_err = vector()
for (i in (1:nrow(rankgrid))) {
  ranktemp = rankgrid[i,]
  try = try(Highdim_AQR(y, X, X_vec, a0, rank = ranktemp, tau = tau, lambda = 0.02,
                        penalty = c(TRUE,TRUE,TRUE),tol = 5e-03, iter_max = 50))
  if("try-error" %in%class(try)){cat("error when rank =", ranktemp);next}
  res = y_val- X_val_vec%*%vec(try$A)
  pre_err[i] = mean(res*(tau - (res<0)))
}

rankhat = rankgrid[which.min(pre_err),]

########estimate  ################
Ahat = Highdim_AQR(y, X, X_vec, a0, rank = rankhat, tau = tau, lambda = 0.02,
                      penalty = c(TRUE,TRUE,TRUE),tol = 4e-03, iter_max = 60)


```
