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
