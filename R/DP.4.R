# unlink("", recursive = TRUE)
#' Releasing 4 Differential Private RKHS smoothing means of a dataset
#'
#' This function create 4 DP RKHS smoothing means from an existing data set
#' with known eigenvalues and eigenvectors
#' @param N real vector 4*1, number of observations
#' @param n real vector 4*1, number of grid points
#' @param e.val.x real vector n*1, eigenvalues
#' @param e.vec.x real valued matrix n*N, eigenvectors
#' @param tau range of the uniform distribution in KL expansion
#' @param phi real number, penalty parameter
#' @param mu real vector n*1, initial mean vector
#' @param m positive integer, number of eigenvectors are going to be used
#' @param e.val.z real valued matrix n*N, eigenvectors of noise
#' @param pow smoothing parameter, e.val.x_{i}=i^{-pow}
#' @param ro range parameter in kernel, real number
#' @return f.tilda: DP RKHS smoothing mean,  n*1 real valued vector
#' @return delta: the coefficient of the noise, real number
#' @return f: RKHS smoothing mean, n*1 real valued vector
DP4=function(alpha=rep(1,4),beta=rep(0.1,4),kernel=c("Exp","M3/2","M5/2","Gau"),
              n=rep(100,4),N=rep(25,4),
              phi=rep(0.01,4),ro=rep(0.2,4),tau=rep(0.4,4),
              case=rep(2,4),pow=rep(4,4),mu=matrix(rep(0,4*100),ncol = 4)){

  # kernal : vactor with 4 elements
  # phi,ro : vactor with 4 elements
  A=list("DP1"=NA,"DP2"=NA,"DP3"=NA,"DP4"=NA,
         "alpha"=alpha,"beta"=beta,"kernel"=kernel,"phi"=phi,
         "ro"=ro,"pow"=pow,"mu"=mu,"n"=n,"N"=N,"tau"=tau,"case"=case,"pow"=pow)
  for(i in 1:4){
    A[[i]]=DP(alpha =alpha[i],beta = beta[i],kernel=kernel[i],phi=phi[i],ro=ro[i],
              n=n[i],N=N[i],tau=tau[i],case=case[i],pow=pow[i],mu=mu[,i])
  }
  return(A)
}

# set.seed(1100)
# F=DP4(alpha=rep(1,4),beta=rep(0.1,4),kernel=c("Exp","M3/2","M5/2","Gau"),
#       n=rep(100,4),N=rep(25,4),
#       phi=rep(0.01,4),ro=rep(0.2,4),tau=rep(0.4,4),
#       case=rep(2,4),pow=rep(4,4),mu=matrix(rep(0.1*sin(2*pi*seq(0,1,length=100)),4),ncol = 4))
