# Functions for generating design matrices and calculating variances
# for staircase designs
#
# Kelsey Grantham (kelsey.grantham@monash.edu)
#

library(ggplot2)
library(reshape2)
library(tidyverse)

# Generate staircase design matrix
SCdesmat <- function(S, K=1, pre=1, post=1) {
  # Inputs:
  #  S - number of treatment sequences/clusters
  #  K - number of times each sequence is repeated
  #  pre - number of pre-switch measurement periods
  #  post - number of post-switch measurement periods
  # Output:
  #  Design matrix
  Xsc <- matrix(data=NA, nrow=S, ncol=(S+pre+post-1))
  for(i in 1:S) {
    Xsc[i,i:(i+pre-1)] <- 0
    Xsc[i,(i+pre):(i+pre+post-1)] <- 1
  }
  Xscreps <- Xsc[sort(rep(1:S, K)), ]
  return(Xscreps)
}

stopifnot(colSums(SCdesmat(3, 2, 1, 1), na.rm=TRUE) == c(0, 2, 2, 2))
stopifnot(colSums(SCdesmat(5, 1, 2, 2), na.rm=TRUE)[2] == 0)
stopifnot(colSums(SCdesmat(4, 2, 1, 2), na.rm=TRUE)[6] == 2)

# Calculate multiple-period CRT treatment effect variance
CRTvartheta <- function(m, Xmat, rho0, r, corrtype, pereff) {
  # Inputs:
  #  m - number of subjects per cluster-period
  #  Xmat - design matrix (period effects and treatment sequences)
  #  rho0 - within-period intracluster correlation
  #  r - cluster autocorrelation
  #  corrtype - within-cluster correlation structure type
  #             (0=block-exchangeable, 1=exponential decay)
  #  pereff - time period effect type
  #           ('cat'=categorical period effects, 'lin'=linear period effects)
  # Output:
  #  Variance of treatment effect estimator
  # Assumptions:
  #  Total variance = 1
  
  totalvar <- 1
  sig2CP <- rho0*totalvar
  sig2E <- totalvar - sig2CP
  sig2 <- sig2E/m
  
  Tp <- ncol(Xmat)
  K <- nrow(Xmat)
  Xvec <- as.vector(t(Xmat))
  
  if(pereff=='cat'){
    stackI <- matrix(rep(diag(1,Tp)), nrow=K*Tp, ncol=Tp, byrow=TRUE)
    Zmat <- cbind(stackI[!is.na(Xvec),], Xvec[!is.na(Xvec)])
  }else if(pereff=='lin'){
    stackT <- matrix(1:Tp, nrow=K*Tp, ncol=1, byrow=TRUE)
    stackT<- cbind(rep(1, nrow(stackT)), stackT)
    Zmat <- cbind(stackT[!is.na(Xvec),], Xvec[!is.na(Xvec)])
  }
  
  # Covariance matrix for one cluster, with decay in correlation over time
  if(corrtype==0){
    # Block-exchangeable structure if corrtype==0
    Vi <-diag(sig2 +(1-r)*sig2CP, Tp) + matrix(data=sig2CP*r, nrow=Tp, ncol=Tp)
  }else if(corrtype==1){
  # Exponential decay structure if corrtype==1
    Vi <- diag(sig2,Tp) + sig2CP*(r^abs(matrix(1:Tp,nrow=Tp, ncol=Tp, byrow=FALSE) -
                                        matrix(1:Tp,nrow=Tp, ncol=Tp, byrow=TRUE)))
  }
  # Covariance matrix for all K clusters
  Vall <- kronecker(diag(1,K), Vi)
  Vall <- Vall[!is.na(Xvec),!is.na(Xvec)]
  
  return(solve((t(Zmat)%*%solve(Vall)%*%Zmat))[ncol(Zmat),ncol(Zmat)])
}

VarSClin <- function(m, S, K, rho0, rhou, r){
  a <- (1 + (m-1)*rho0)/m
  b <- (rhou + (m-1)*r*rho0)/m
  
  vartheta <- (2*((S^2+2)*a - (S^2-4)*b))/(K*S*(S^2-1))
  return(vartheta)
}

VarSCcat <- function(m, S, K, rho0, rhou, r){
  a <- (1 + (m-1)*rho0)/m
  b <- (rhou + (m-1)*r*rho0)/m

  fracterm <- ((a + sqrt(a^2-b^2))^S - b^S)/((a + sqrt(a^2-b^2))^S + b^S)
  vartheta <- (2*(a-b)^2)/(K*(S*(a-b)-sqrt(a^2-b^2)*fracterm))
  return(vartheta)
}

pow <- function(vars, effsize, siglevel=0.05){
  z <- qnorm(siglevel/2)
  pow <- pnorm(z + sqrt(1/vars)*effsize)
  return(pow)
}
