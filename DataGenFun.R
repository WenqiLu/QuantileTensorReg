########################
### Generate coef A
#######################
DGP_A = function(sigma_g, dim = c(10,10,5), rank = c(3,3,3)){
  ###nonsparse
  ###sigma_g: sd of core tensor G
  ###dim: dimension of tensor A
  ###rank: rank of tensor A
  ###U: list of factor matrix
  ###G: core tensor
  order = length(dim)
  U = list()
  for (i in 1:order) {
    Ui = matrix(rnorm(dim[i]*dim[i]),dim[i],dim[i])
    Ui = Ui%*%t(Ui)
    U[[i]] = eigen(Ui)$vectors[,1:rank[i]]
  }
  G = rand_tensor(modes = rank)
  A = ttm(G,U[[1]],m = 1)
  for (i in 2:order) {
    A = ttm(A,U[[i]],m = i)
  }
  rtn = list(A = A, U = U)
  return(rtn)
}


DGP_Ahigh = function(sigma_g, dim = c(10,10,10), rank = c(2,2,2)){
  ###sparse A
  ###sigma_g: sd of core tensor G
  ###dim: dimension of tensor A
  ###rank: rank of tensor A
  ###U: list of factor matrix
  ###G: core tensor
  order = length(dim)
  U = list()
  
  i=1
  nn=3
  U[[i]] = cbind(rnorm(nn),rep(0,nn))
  U[[i]] = rbind(U[[i]],cbind(rep(0,5-nn),rnorm(5-nn)))
  for(ii in 2:(dim[i]/5)){
    U[[i]] = rbind(U[[i]],cbind(rnorm(nn),rep(0,nn)))
    U[[i]] = rbind(U[[i]],cbind(rep(0,5-nn),rnorm(5-nn)))
  }
  U[[i]] = apply(U[[i]], 2, function(x) x/sqrt(sum(x^2)))
  if(all(rank ==c(2,2,2))){
    for(i in 2:3){
      nn=3
      U[[i]] = cbind(rnorm(nn),rep(0,nn))
      U[[i]] = rbind(U[[i]],cbind(rep(0,5-nn),rnorm(5-nn)))
      for(ii in 2:(dim[i]/5)){
        U[[i]] = rbind(U[[i]],cbind(rnorm(nn),rep(0,nn)))
        U[[i]] = rbind(U[[i]],cbind(rep(0,5-nn),rnorm(5-nn)))
      }
      U[[i]] = apply(U[[i]], 2, function(x) x/sqrt(sum(x^2)))
    }
  }else{
    uu = matrix(rnorm(nn^2),nn,nn)
    uu = eigen(t(uu)%*%uu)$vectors[,1:2]
    U[[2]] = cbind(uu,rep(0,nn))
    bb = rnorm(5-nn)
    bb = bb/sqrt(sum(bb^2))
    U[[2]] = rbind(U[[2]],cbind(matrix(0,5-nn,2),bb))
    for (ii in 2:(dim[2]/5)) {
      uu = matrix(rnorm(nn^2),nn,nn)
      uu = eigen(t(uu)%*%uu)$vectors[,1:2]
      U[[2]] = rbind(U[[2]],cbind(uu,rep(0,nn)))  
      bb = rnorm(5-nn)
      bb = bb/sqrt(sum(bb^2))
      U[[2]] = rbind(U[[2]],cbind(matrix(0,5-nn,2),bb))
    }
    U[[2]] = apply(U[[2]], 2, function(x) x/sqrt(sum(x^2)))
    
    uu = matrix(0,nn,2)
    bb = matrix(rnorm(nn^2),nn,nn)
    bb = eigen(t(bb)%*%bb)$vectors[,1:2]
    U[[3]] = cbind(uu,bb)
    uu = matrix(rnorm((5-nn)^2),5-nn,5-nn)
    uu = eigen(t(uu)%*%uu)$vectors[,1:2]
    bb = matrix(0,5-nn,2)
    U[[3]] = rbind(U[[3]],cbind(uu,bb))  
    for (ii in 2:(dim[3]/5)) {
      uu = matrix(0,nn,2)
      bb = matrix(rnorm(nn^2),nn,nn)
      bb = eigen(t(bb)%*%bb)$vectors[,1:2]
      U[[3]] = rbind(U[[3]],cbind(uu,bb))
      uu = matrix(rnorm((5-nn)^2),5-nn,5-nn)
      uu = eigen(t(uu)%*%uu)$vectors[,1:2]
      bb = matrix(0,5-nn,2)
      U[[3]] = rbind(U[[3]],cbind(uu,bb))  
    }
    U[[3]] = apply(U[[3]], 2, function(x) x/sqrt(sum(x^2)))
  }
   
  G = rand_tensor(modes = rank)
  A = ttm(G,U[[1]],m = 1)
  for (i in 2:order) {
    A = ttm(A,U[[i]],m = i)
  }
  rtn = list(A = A, U = U, G = G)
  return(rtn)
}

######################
### Generate X, y
######################
DGP = function(n, coef, sigma, tau, dim = c(10,10,10)) {
  ###model: Y_i=<A,X_i>+ep_i
  ###sigma: sd of error term epsilon (ep)
  ###tau: tau th quantile
  ###dim: dimension of tensor A
  ###rank: rank of tensor A
  signal = vector(length = n)
  X = list()
  for (i in 1:n) {
    X_i = rand_tensor(modes = dim)
    X[[i]] = X_i
    signal[i] = vec(coef)%*%vec(X_i)
    #print(i)
  }
  q_tau = qnorm(tau,0,1)
  ep = rnorm(n,mean = -q_tau, sd = sigma)
  y = signal + ep
  rtn = list(X = X, y = y)
  return(rtn)
}

X_unfold = function(tnsrls, mode){
  #vectorize tensor X_i by vectorize mode-k matricization of X_i---vec(X_(k))
  #and storage in a n-row matrix
  ###tnsrls: a list of tensor
  n = length(tnsrls)
  dim = tnsrls[[1]]@modes
  xvec = k_unfold(tnsrls[[1]],mode)
  X_vec = matrix(xvec@data,nrow = 1)
  for (i in 2:n) {
    xvec = k_unfold(tnsrls[[i]],mode)
    X_vec = rbind(X_vec,matrix(xvec@data,nrow = 1))
  }
  return(X_vec)
}