# Efficiency comparisons for staircase design
#
# Kelsey Grantham (kelsey.grantham@monash.edu)
#
# Some functions based on previous work by J Kasza

library(ggplot2)
library(reshape2)

# Generate stepped wedge design matrix
# with K sequences and reps repeated sequences
SWdesmat <- function(K, reps=1) {
  # Inputs:
  #  K - number of unique treatment sequences
  #  reps - number of times each sequence is repeated
  #  (K*reps = number of clusters)
  # Output:
  #  Design matrix
  Xsw <- matrix(data=0, ncol=(K+1), nrow=K)
  for(i in 1:K) {
    Xsw[i,(i+1):(K+1)] <- 1
  }
  return(Xsw[sort(rep(1:K, reps)), ])
}

# Generate staircase design matrix
SCdesmatprepost <- function(K, pre=1, post=1) {
  # Inputs:
  #  K - number of treatment sequences/clusters
  #  pre - number of pre-switch measurement periods
  #  post - number of post-switch measurement periods
  # Output:
  #  Design matrix
  Xsc <- matrix(data=NA, nrow=K, ncol=(K+pre+post-1))
  for(i in 1:K) {
    Xsc[i,i:(i + pre-1)] <- 0
    Xsc[i,(i+pre):(i+pre+post-1)] <- 1
  }
  return(Xsc)
}

# Calculate multiple-period CRT treatment effect variance
CRTVarSW <- function(m, Xmat, rho0, r, corrtype, pereff) {
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

varSWSC_plot <- function(m, S, T0, T1, corrtype, pereff){
  # Compare variances of complete SW and staircase designs, each with S
  # sequences/clusters, for a range of correlation parameters
  # Inputs:
  #  m - number of subjects per cluster-period
  #  S - number of treatment sequences/clusters
  #  T0 - number of pre-switch measurement periods
  #  T1 - number of post-switch measurement periods
  #  corrtype - within-cluster correlation structure type
  #             (0=block-exchangeable, 1=exponential decay)
  #  pereff - time period effect type
  #           ('cat'=categorical period effects, 'lin'=linear period effects)
  # Output:
  #  Contour plot of relative variances (vartheta_SC/vartheta_SW)
  
  rho0seq <- seq(0.01, 0.99, 0.05)
  rseq <- seq(0.0, 0.95, 0.05)
  
  SWSCvars <- matrix(data=NA, nrow= length(rho0seq), ncol=length(rseq))
  for(i in 1:length(rho0seq)) {
    for(rind in 1:length(rseq)) {
      SWSCvars[i,rind] <-CRTVarSW(m, SCdesmatprepost(S, T0, T1), rho0seq[i],
                                  rseq[rind], corrtype=corrtype, pereff=pereff)/
                          CRTVarSW(m, SWdesmat(S), rho0seq[i],
                                   rseq[rind], corrtype=corrtype, pereff=pereff)
    }
  }
  
  # Plot the results using a contour plot
  SWSCvars<-round(SWSCvars, 2)
  meltSWSCvars <- melt(SWSCvars)
  
  names(meltSWSCvars)[names(meltSWSCvars)=="Var1"] <- "rho"
  names(meltSWSCvars)[names(meltSWSCvars)=="Var2"] <- "r"
  
  rhovec <- as.vector(matrix(data=rho0seq, nrow=length(rho0seq), ncol=length(rseq), byrow=FALSE))
  rvec <- as.vector(matrix(data=rseq, nrow=length(rho0seq), ncol=length(rseq), byrow=TRUE))
  meltSWSCvars$rhoseq <- rhovec
  meltSWSCvars$rseq <- rvec

  myplot <- ggplot(meltSWSCvars, aes(x=rseq, y=rhoseq)) + 
    geom_tile(aes(fill=value)) + 
    scale_fill_gradientn(colours=c("yellow","red")) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    theme(aspect.ratio=3/8, legend.position="none", legend.key.size=unit(1, "cm"), 
          legend.text=element_text(size=12), 
          legend.background = element_rect(fill="grey95")) +
    coord_fixed() + xlab("Cluster autocorrelation, r") +  ylab("Within-period ICC") +
    geom_text(aes(rseq, rhoseq, label=round(value,2)), color="black", size=3) 

  return(myplot)
}
