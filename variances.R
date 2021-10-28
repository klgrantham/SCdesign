# Functions for generating design matrices, calculating variances and plotting
# relative variances, for complete stepped wedge and staircase designs
#
# Kelsey Grantham (kelsey.grantham@monash.edu)
#
# Some functions based on previous work by J Kasza

library(ggplot2)
library(reshape2)

# Generate stepped wedge design matrix
# with K sequences and reps repeated sequences
SWdesmat <- function(S, reps=1) {
  # Inputs:
  #  S - number of unique treatment sequences
  #  reps - number of times each sequence is repeated
  #  (S*reps = number of clusters)
  # Output:
  #  Design matrix
  Xsw <- matrix(data=0, ncol=(S+1), nrow=S)
  for(i in 1:S) {
    Xsw[i,(i+1):(S+1)] <- 1
  }
  Xswreps <- Xsw[sort(rep(1:S, reps)), ]
  return(Xswreps)
}

stopifnot(colSums(SWdesmat(3, 1))[4] == 3)
stopifnot(colSums(SWdesmat(5, 2))[3] == 4)

# Generate staircase design matrix
SCdesmat <- function(S, pre=1, post=1, reps=1) {
  # Inputs:
  #  S - number of treatment sequences/clusters
  #  reps - number of times each sequence is repeated
  #  pre - number of pre-switch measurement periods
  #  post - number of post-switch measurement periods
  # Output:
  #  Design matrix
  Xsc <- matrix(data=NA, nrow=S, ncol=(S+pre+post-1))
  for(i in 1:S) {
    Xsc[i,i:(i+pre-1)] <- 0
    Xsc[i,(i+pre):(i+pre+post-1)] <- 1
  }
  Xscreps <- Xsc[sort(rep(1:S, reps)), ]
  return(Xscreps)
}

stopifnot(colSums(SCdesmat(3, 1, 1, 2), na.rm=TRUE) == c(0, 2, 2, 2))
stopifnot(colSums(SCdesmat(5, 2, 2, 1), na.rm=TRUE)[2] == 0)
stopifnot(colSums(SCdesmat(4, 1, 2, 2), na.rm=TRUE)[6] == 2)

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

varSCSW_plot_general <- function(m_SW, S_SW, reps_SW, m_SC, S_SC, pre_SC, post_SC,
                                 reps_SC, corrtype, pereff){
  # Compare variances of complete SW and staircase designs, for a range of
  # correlation parameters
  # Inputs:
  #  m_SW - number of subjects per cluster-period for SW design
  #  S_SW - number of unique treatment sequences for SW design
  #  reps_SW - number of times each sequence is repeated for SW design
  #  m_SC - number of subjects per cluster-period for SC design
  #  S_SC - number of unique treatment sequences for SC design
  #  pre_SC - number of pre-switch measurement periods for SC design
  #  post_SC - number of post-switch measurement periods for SC design
  #  corrtype - within-cluster correlation structure type
  #             (0=block-exchangeable, 1=exponential decay)
  #  pereff - time period effect type
  #           ('cat'=categorical period effects, 'lin'=linear period effects)
  # Output:
  #  Contour plot of relative variances (vartheta_SC/vartheta_SW)
  
  rho0seq <- seq(0.01, 0.99, 0.05)
  rseq <- seq(0.0, 0.95, 0.05)
  
  SWSCvars <- matrix(data=NA, nrow=length(rho0seq), ncol=length(rseq))
  for(i in 1:length(rho0seq)) {
    for(rind in 1:length(rseq)) {
      SWSCvars[i,rind] <-CRTVarSW(m_SC, SCdesmat(S_SC, pre_SC, post_SC, reps_SC),
                      rho0seq[i], rseq[rind], corrtype=corrtype, pereff=pereff)/
                        CRTVarSW(m_SW, SWdesmat(S_SW, reps_SW), rho0seq[i],
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
