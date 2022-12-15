# Functions for generating design matrices, calculating variances and plotting
# relative variances, for complete stepped wedge and staircase designs
#
# Kelsey Grantham (kelsey.grantham@monash.edu)
#
# Some functions based on previous work by J Kasza

library(ggplot2)
library(reshape2)
library(tidyverse)
library(viridis)

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
SCdesmat <- function(S, reps=1, pre=1, post=1) {
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

stopifnot(colSums(SCdesmat(3, 2, 1, 1), na.rm=TRUE) == c(0, 2, 2, 2))
stopifnot(colSums(SCdesmat(5, 1, 2, 2), na.rm=TRUE)[2] == 0)
stopifnot(colSums(SCdesmat(4, 2, 1, 2), na.rm=TRUE)[6] == 2)

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

gridvals_small <- function(m_SW, S_SW, reps_SW, m_SC, S_SC, reps_SC,
                     pre_SC, post_SC, corrtype, pereff){
  # Compare variances of complete SW and staircase designs, for a range of
  # correlation parameters
  # Inputs:
  #  m_SW - number of subjects per cluster-period for SW design
  #  S_SW - number of unique treatment sequences for SW design
  #  reps_SW - number of times each sequence is repeated for SW design
  #  m_SC - number of subjects per cluster-period for SC design
  #  S_SC - number of unique treatment sequences for SC design
  #  reps_SC - number of times each sequence is repeated for SC design
  #  pre_SC - number of pre-switch measurement periods for SC design
  #  post_SC - number of post-switch measurement periods for SC design
  #  corrtype - within-cluster correlation structure type
  #             (0=block-exchangeable, 1=exponential decay)
  #  pereff - time period effect type
  #           ('cat'=categorical period effects, 'lin'=linear period effects)
  # Output:
  #  Relative variances (vartheta_SC/vartheta_SW) for a range of rho and r values
  
  rho0seq <- seq(0.01, 0.25, 0.01)
  rseq <- seq(0.2, 0.95, 0.05)
  
  SCSWvars <- matrix(data=NA, nrow=length(rho0seq), ncol=length(rseq))
  for(i in 1:length(rho0seq)) {
    for(rind in 1:length(rseq)) {
      SCSWvars[i,rind] <-CRTVarSW(m_SC, SCdesmat(S_SC, reps_SC, pre_SC, post_SC),
                                  rho0seq[i], rseq[rind], corrtype=corrtype, pereff=pereff)/
        CRTVarSW(m_SW, SWdesmat(S_SW, reps_SW), rho0seq[i],
                 rseq[rind], corrtype=corrtype, pereff=pereff)
    }
  }
  
  # Plot the results using a contour plot
  SCSWvars<-round(SCSWvars, 2)
  meltSCSWvars <- melt(SCSWvars)
  
  names(meltSCSWvars)[names(meltSCSWvars)=="Var1"] <- "rho"
  names(meltSCSWvars)[names(meltSCSWvars)=="Var2"] <- "r"
  
  rhovec <- as.vector(matrix(data=rho0seq, nrow=length(rho0seq), ncol=length(rseq), byrow=FALSE))
  rvec <- as.vector(matrix(data=rseq, nrow=length(rho0seq), ncol=length(rseq), byrow=TRUE))
  meltSCSWvars$rhoseq <- rhovec
  meltSCSWvars$rseq <- rvec
  return(meltSCSWvars)
}

varSCSW_grid_small_plot <- function(m_SW, S_SW, reps_SW, m_SC, S_SC, reps_SC,
                              pre_SC, post_SC, corrtype, pereff){
  # Compare variances of complete SW and staircase designs, for a range of
  # correlation parameters
  # Inputs:
  #  m_SW - number of subjects per cluster-period for SW design
  #  S_SW - number of unique treatment sequences for SW design
  #  reps_SW - number of times each sequence is repeated for SW design
  #  m_SC - number of subjects per cluster-period for SC design
  #  S_SC - number of unique treatment sequences for SC design
  #  reps_SC - number of times each sequence is repeated for SC design
  #  pre_SC - number of pre-switch measurement periods for SC design
  #  post_SC - number of post-switch measurement periods for SC design
  #  corrtype - within-cluster correlation structure type
  #             (0=block-exchangeable, 1=exponential decay)
  #  pereff - time period effect type
  #           ('cat'=categorical period effects, 'lin'=linear period effects)
  # Output:
  #  Contour plot of relative variances (vartheta_SC/vartheta_SW)
  
  relvars <- gridvals_small(m_SW, S_SW, reps_SW, m_SC, S_SC, reps_SC, pre_SC, post_SC,
                      corrtype, pereff)  
  
  p <- ggplot(relvars, aes(x=rseq, y=rhoseq)) + 
    geom_tile(aes(fill=value)) + 
    scale_fill_gradientn(colours=c("yellow","red")) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    theme(aspect.ratio=3/8, legend.position="bottom", legend.key.size=unit(1, "cm"), 
          legend.text=element_text(size=12), 
          legend.background = element_rect(fill="grey95"),
          plot.title=element_text(hjust=0.5, size=16),
          axis.title=element_text(size=14), axis.text=element_text(size=14)) +
    coord_fixed() + xlab(expression(paste("Cluster autocorrelation, ", r))) + ylab(expression(paste("Within-period ICC, ", rho))) +
    geom_text(aes(rseq, rhoseq, label=round(value,2)), color="black", size=3)  +
    ggtitle(bquote(paste("Relative variance of treatment effect estimators, ",
                        Var(hat(theta))[paste("SC(", .(S_SC), ",", .(reps_SC), ",", .(pre_SC), ",", .(post_SC), ")")]/
                        Var(hat(theta))[paste("SW(", .(S_SW), ",", .(reps_SW), ")")])))
  corrname <- ifelse(corrtype==0, "BE", "DTD")
  ggsave(paste0("plots/SC_", S_SC, reps_SC, pre_SC, post_SC, "_vs_SW_",
                S_SW, reps_SW, "_", corrname, "_", pereff, ".jpg"),
         p, width=9, height=7, units="in", dpi=800)
  return(p)
}

varSCSW_grid_multiplot <- function(S_SW, reps_SW, S_SC, reps_SC,
                                    pre_SC, post_SC, corrtype){
  # Compare variances of complete SW and staircase designs, for a range of
  # correlation parameters
  # Inputs:
  #  m_SW - number of subjects per cluster-period for SW design
  #  S_SW - number of unique treatment sequences for SW design
  #  reps_SW - number of times each sequence is repeated for SW design
  #  m_SC - number of subjects per cluster-period for SC design
  #  S_SC - number of unique treatment sequences for SC design
  #  reps_SC - number of times each sequence is repeated for SC design
  #  pre_SC - number of pre-switch measurement periods for SC design
  #  post_SC - number of post-switch measurement periods for SC design
  #  corrtype - within-cluster correlation structure type
  #             (0=block-exchangeable, 1=exponential decay)
  #  pereff - time period effect type
  #           ('cat'=categorical period effects, 'lin'=linear period effects)
  # Output:
  #  Contour plot of relative variances (vartheta_SC/vartheta_SW)
  
  m1 <- 10
  m2 <- 100
  relvars_m1_cat <- gridvals_small(m1, S_SW, reps_SW, m1, S_SC, reps_SC,
                                   pre_SC, post_SC, corrtype, 'cat')
  relvars_m1_cat$m <- "m = 10"
  relvars_m1_cat$pereff <- "Categorical time"
  relvars_m2_cat <- gridvals_small(m2, S_SW, reps_SW, m2, S_SC, reps_SC,
                                   pre_SC, post_SC, corrtype, 'cat')
  relvars_m2_cat$m <- "m = 100"
  relvars_m2_cat$pereff <- "Categorical time"
  relvars_m1_lin <- gridvals_small(m1, S_SW, reps_SW, m1, S_SC, reps_SC,
                                   pre_SC, post_SC, corrtype, 'lin')
  relvars_m1_lin$m <- "m = 10"
  relvars_m1_lin$pereff <- "Linear time"
  relvars_m2_lin <- gridvals_small(m2, S_SW, reps_SW, m2, S_SC, reps_SC,
                                   pre_SC, post_SC, corrtype, 'lin')
  relvars_m2_lin$m <- "m = 100"
  relvars_m2_lin$pereff <- "Linear time"
  relvars <- bind_rows(
    relvars_m1_cat,
    relvars_m2_cat,
    relvars_m1_lin,
    relvars_m2_lin
  )
  
  p <- ggplot(relvars, aes(x=rseq, y=rhoseq)) + 
    geom_tile(aes(fill=value)) +
    facet_grid(
      m ~ pereff
    ) +
    scale_fill_viridis_c(name="Relative variance", direction=-1) +
#    scale_fill_gradientn(colours=c("yellow","red")) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    theme(aspect.ratio=3/8,
          legend.key.width = unit(1.5, "cm"),
          legend.title=element_text(size=14), legend.text=element_text(size=14),
          legend.position="bottom",
          plot.title=element_text(hjust=0.5, size=16),
          axis.title=element_text(size=14), axis.text=element_text(size=14),
          strip.background = element_rect(
            color="white", fill="white", linetype="solid"
          ),
          strip.text.x = element_text(size = 12),
          strip.text.y = element_text(size=12)) +
    coord_fixed() + xlab(expression(paste("Cluster autocorrelation, ", r))) + ylab(expression(paste("Within-period ICC, ", rho))) +
#    geom_text(aes(rseq, rhoseq, label=round(value,2)), color="black", size=3)  +
    ggtitle(bquote(paste("Relative variance of treatment effect estimators, ",
                         Var(hat(theta))[paste("SC(", .(S_SC), ",", .(reps_SC), ",", .(pre_SC), ",", .(post_SC), ")")]/
                           Var(hat(theta))[paste("SW(", .(S_SW), ",", .(reps_SW), ")")])))
  corrname <- ifelse(corrtype==0, "BE", "DTD")
  ggsave(paste0("plots/multiplot_SC_", S_SC, reps_SC, pre_SC, post_SC, "_vs_SW_",
                S_SW, reps_SW, "_", corrname, ".jpg"),
         p, width=9, height=5, units="in", dpi=800)
  return(p)
}

gridvals_varSC <- function(m_SC, S_SC, reps_SC, pre_SC, post_SC, corrtype, pereff){
  # Get variances of staircase design, for a range of correlation parameters
  # Inputs:
  #  m_SC - number of subjects per cluster-period for SC design
  #  S_SC - number of unique treatment sequences for SC design
  #  reps_SC - number of times each sequence is repeated for SC design
  #  pre_SC - number of pre-switch measurement periods for SC design
  #  post_SC - number of post-switch measurement periods for SC design
  #  corrtype - within-cluster correlation structure type
  #             (0=block-exchangeable, 1=exponential decay)
  #  pereff - time period effect type
  #           ('cat'=categorical period effects, 'lin'=linear period effects)
  # Output:
  #  Variances for a range of rho and r values
  
  rho0seq <- seq(0.05, 0.95, 0.05)
  rseq <- seq(0.05, 0.95, 0.05)
  
  SCvars <- matrix(data=NA, nrow=length(rho0seq), ncol=length(rseq))
  for(i in 1:length(rho0seq)) {
    for(rind in 1:length(rseq)) {
      SCvars[i,rind] <-CRTVarSW(m_SC, SCdesmat(S_SC, reps_SC, pre_SC, post_SC),
                                rho0seq[i], rseq[rind], corrtype=corrtype, pereff=pereff)
    }
  }
  
  # Plot the results using a contour plot
  SCvars<-round(SCvars, 4)
  meltSCvars <- melt(SCvars)
  
  names(meltSCvars)[names(meltSCvars)=="Var1"] <- "rho"
  names(meltSCvars)[names(meltSCvars)=="Var2"] <- "r"
  
  rhovec <- as.vector(matrix(data=rho0seq, nrow=length(rho0seq), ncol=length(rseq), byrow=FALSE))
  rvec <- as.vector(matrix(data=rseq, nrow=length(rho0seq), ncol=length(rseq), byrow=TRUE))
  meltSCvars$rhoseq <- rhovec
  meltSCvars$rseq <- rvec
  return(meltSCvars)
}

varSC_grid_multiplot <- function(S1_SC, S2_SC, reps_SC, pre_SC, post_SC, corrtype, pereff){
  # Display variances of staircase designs, for a range of
  # correlation parameters
  # Inputs:
  #  S1_SC - number of unique treatment sequences for SC design (left column)
  #  S2_SC - number of unique treatment sequences for SC design (right column)
  #  reps_SC - number of times each sequence is repeated for SC design
  #  pre_SC - number of pre-switch measurement periods for SC design
  #  post_SC - number of post-switch measurement periods for SC design
  #  corrtype - within-cluster correlation structure type
  #             (0=block-exchangeable, 1=exponential decay)
  #  pereff - time period effect type
  #           ('cat'=categorical period effects, 'lin'=linear period effects)
  # Output:
  #  Contour plot of relative variances (vartheta_SC/vartheta_SW)
  
  m1 <- 10
  m2 <- 100
  S1 <- paste("S = ", S1_SC)
  S2 <- paste("S = ", S2_SC)
  
  vars_m1_S1 <- gridvals_varSC(m1, S1_SC, reps_SC, pre_SC, post_SC, corrtype, pereff)
  vars_m1_S1$m <- "m = 10"
  vars_m1_S1$Sname <- S1
  vars_m2_S1 <- gridvals_varSC(m2, S1_SC, reps_SC, pre_SC, post_SC, corrtype, pereff)
  vars_m2_S1$m <- "m = 100"
  vars_m2_S1$Sname <- S1
  vars_m1_S2 <- gridvals_varSC(m1, S2_SC, reps_SC, pre_SC, post_SC, corrtype, pereff)
  vars_m1_S2$m <- "m = 10"
  vars_m1_S2$Sname <- S2
  vars_m2_S2 <- gridvals_varSC(m2, S2_SC, reps_SC, pre_SC, post_SC, corrtype, pereff)
  vars_m2_S2$m <- "m = 100"
  vars_m2_S2$Sname <- S2
  vars <- bind_rows(
    vars_m1_S1,
    vars_m2_S1,
    vars_m1_S2,
    vars_m2_S2
  )
  
  p <- ggplot(vars, aes(x=rseq, y=rhoseq)) + 
    geom_tile(aes(fill=value)) +
    facet_grid(
      m ~ Sname
    ) +
    scale_fill_viridis_c(name="Variance", direction=-1) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    theme(aspect.ratio=3/8,
          legend.key.width = unit(1.5, "cm"),
          legend.title=element_text(size=14), legend.text=element_text(size=14),
          legend.position="bottom",
          plot.title=element_text(hjust=0.5, size=16),
          axis.title=element_text(size=14), axis.text=element_text(size=14),
          strip.background = element_rect(
            color="white", fill="white", linetype="solid"
          ),
          strip.text.x = element_text(size = 12),
          strip.text.y = element_text(size=12)) +
    coord_fixed() + xlab(expression(paste("Cluster autocorrelation, ", r))) + ylab(expression(paste("Within-period ICC, ", rho))) +
    ggtitle(bquote(paste("Variance of treatment effect estimator, ",
                         Var(hat(theta))[paste("SC(S,", .(reps_SC), ",", .(pre_SC), ",", .(post_SC), ")")])))
  ggsave(paste0("plots/multiplot_SC_", S1_SC, "vs", S2_SC, reps_SC, pre_SC, post_SC, "_", corrtype, "_", pereff, ".jpg"),
         p, width=9, height=5, units="in", dpi=800)
  return(p)
}

gridvals <- function(m_SW, S_SW, reps_SW, m_SC, S_SC, reps_SC,
                     pre_SC, post_SC, corrtype, pereff){
  # Compare variances of complete SW and staircase designs, for a range of
  # correlation parameters
  # Inputs:
  #  m_SW - number of subjects per cluster-period for SW design
  #  S_SW - number of unique treatment sequences for SW design
  #  reps_SW - number of times each sequence is repeated for SW design
  #  m_SC - number of subjects per cluster-period for SC design
  #  S_SC - number of unique treatment sequences for SC design
  #  reps_SC - number of times each sequence is repeated for SC design
  #  pre_SC - number of pre-switch measurement periods for SC design
  #  post_SC - number of post-switch measurement periods for SC design
  #  corrtype - within-cluster correlation structure type
  #             (0=block-exchangeable, 1=exponential decay)
  #  pereff - time period effect type
  #           ('cat'=categorical period effects, 'lin'=linear period effects)
  # Output:
  #  Relative variances (vartheta_SC/vartheta_SW) for a range of rho and r values
  
  rho0seq <- seq(0.01, 0.99, 0.05)
  rseq <- seq(0.0, 0.95, 0.05)
  
  SCSWvars <- matrix(data=NA, nrow=length(rho0seq), ncol=length(rseq))
  for(i in 1:length(rho0seq)) {
    for(rind in 1:length(rseq)) {
      SCSWvars[i,rind] <-CRTVarSW(m_SC, SCdesmat(S_SC, reps_SC, pre_SC, post_SC),
                      rho0seq[i], rseq[rind], corrtype=corrtype, pereff=pereff)/
                         CRTVarSW(m_SW, SWdesmat(S_SW, reps_SW), rho0seq[i],
                      rseq[rind], corrtype=corrtype, pereff=pereff)
    }
  }
  
  # Plot the results using a contour plot
  SCSWvars<-round(SCSWvars, 2)
  meltSCSWvars <- melt(SCSWvars)
  
  names(meltSCSWvars)[names(meltSCSWvars)=="Var1"] <- "rho"
  names(meltSCSWvars)[names(meltSCSWvars)=="Var2"] <- "r"
  
  rhovec <- as.vector(matrix(data=rho0seq, nrow=length(rho0seq), ncol=length(rseq), byrow=FALSE))
  rvec <- as.vector(matrix(data=rseq, nrow=length(rho0seq), ncol=length(rseq), byrow=TRUE))
  meltSCSWvars$rhoseq <- rhovec
  meltSCSWvars$rseq <- rvec
  return(meltSCSWvars)
}

varSCSW_grid_plot <- function(m_SW, S_SW, reps_SW, m_SC, S_SC, reps_SC,
                              pre_SC, post_SC, corrtype, pereff){
  # Compare variances of complete SW and staircase designs, for a range of
  # correlation parameters
  # Inputs:
  #  m_SW - number of subjects per cluster-period for SW design
  #  S_SW - number of unique treatment sequences for SW design
  #  reps_SW - number of times each sequence is repeated for SW design
  #  m_SC - number of subjects per cluster-period for SC design
  #  S_SC - number of unique treatment sequences for SC design
  #  reps_SC - number of times each sequence is repeated for SC design
  #  pre_SC - number of pre-switch measurement periods for SC design
  #  post_SC - number of post-switch measurement periods for SC design
  #  corrtype - within-cluster correlation structure type
  #             (0=block-exchangeable, 1=exponential decay)
  #  pereff - time period effect type
  #           ('cat'=categorical period effects, 'lin'=linear period effects)
  # Output:
  #  Contour plot of relative variances (vartheta_SC/vartheta_SW)

  relvars <- gridvals(m_SW, S_SW, reps_SW, m_SC, S_SC, reps_SC, pre_SC, post_SC,
                      corrtype, pereff)  

  myplot <- ggplot(relvars, aes(x=rseq, y=rhoseq)) + 
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

varSCSW_line_plot <- function(m_SW, S_SW, reps_SW, m_SC, S_SC, reps_SC,
                              pre_SC, post_SC, corrtype, pereff, title=""){
  
  rhovals <- seq(0.01, 0.2, 0.005)
  rvals <- c(0.25, 0.5, 0.75, 0.95, 1.0)
  relvars <- expand.grid(rho=rhovals, r=rvals)
  relvars$varSC <- with(
    relvars,
    sapply(1:nrow(relvars), function(j){
      CRTVarSW(m_SC, SCdesmat(S_SC, reps_SC, pre_SC, post_SC),
               rho[j], r[j], corrtype, pereff)
      }
    )
  )
  relvars$varSW <- with(
    relvars,
    sapply(1:nrow(relvars), function(j){
      CRTVarSW(m_SW, SWdesmat(S_SW, reps_SW),
               rho[j], r[j], corrtype, pereff)
      }
    )
  )
  relvars <- relvars %>%
    mutate(relvarSCSW = varSC/varSW,
           rfac=as.factor(r)) %>%
    select(c(rho, rfac, relvarSCSW))

  p <- ggplot(data=relvars, aes(x=rho, y=relvarSCSW, colour=rfac)) +
    geom_line(size=1.2) +
    geom_hline(yintercept=1, linetype="dashed") +
#      expand_limits(y=ylimits) +
    xlab("Within-period ICC") +
    ylab("Relative variance") +
    labs(title=title, colour="Cluster autocorrelation") +
    theme_bw() +
    theme(plot.title=element_text(hjust=0.5, size=12),
          axis.title=element_text(size=10), axis.text=element_text(size=10),
          legend.key.width = unit(1.5, "cm"),
          legend.title=element_text(size=12), legend.text=element_text(size=12),
          legend.position="bottom")
  return(p)
}

varSCSW_multi_plot <- function(S_SW, reps_SW, S_SC, reps_SC, pre_SC, post_SC, corrtype){
  # Compare variances of complete SW and staircase designs, for a range of
  # correlation parameters
  # Inputs:
  #  S_SW - number of unique treatment sequences for SW design
  #  reps_SW - number of times each sequence is repeated for SW design
  #  S_SC - number of unique treatment sequences for SC design
  #  reps_SC - number of times each sequence is repeated for SC design
  #  pre_SC - number of pre-switch measurement periods for SC design
  #  post_SC - number of post-switch measurement periods for SC design
  #  corrtype - within-cluster correlation structure type
  #             (0=block-exchangeable, 1=exponential decay)
  # Output:
  #  Multiplot of relative variances (vartheta_SC/vartheta_SW), for m=10 and 100
  #  and categorical and linear period effects
  
  m1 <- 10
  m2 <- 100
  
  # m=10, categorical period effects
  p1 <- varSCSW_line_plot(
    m1, S_SW, reps_SW,
    m1, S_SC, reps_SC, pre_SC, post_SC,
    corrtype, 'cat',
    bquote(paste("m = ", .(m1), ", ", "categorical period effects"))
    )
  # m=100, categorical period effects
  p2 <- varSCSW_line_plot(
    m2, S_SW, reps_SW,
    m2, S_SC, reps_SC, pre_SC, post_SC,
    corrtype, 'cat',
    bquote(paste("m = ", .(m2), ", ", "categorical period effects"))
    )
  # m=10, linear period effects
  p3 <- varSCSW_line_plot(
    m1, S_SW, reps_SW,
    m1, S_SC, reps_SC, pre_SC, post_SC,
    corrtype, 'lin',
    bquote(paste("m = ", .(m1), ", ", "linear period effects"))
    )
  # m=100, linear period effects
  p4 <- varSCSW_line_plot(
    m2, S_SW, reps_SW,
    m2, S_SC, reps_SC, pre_SC, post_SC,
    corrtype, 'lin',
    bquote(paste("m = ", .(m2), ", ", "linear period effects"))
    )
  mylegend <- g_legend(p1)
  title <- bquote(paste("Relative variance of treatment effect estimators, ",
                        Var(hat(theta))[paste("SC(", .(S_SC), ",", .(reps_SC), ",", .(pre_SC), ",", .(post_SC), ")")]/
                        Var(hat(theta))[paste("SW(", .(S_SW), ",", .(reps_SW), ")")]))
  p1to4 <- make_2x2_multiplot(p1, p2, p3, p4, mylegend, title=title)
  corrname <- ifelse(corrtype==0, "BE", "DTD")
  ggsave(paste0("plots/SC", S_SC, reps_SC, pre_SC, post_SC, "_vs_SW",
                S_SW, reps_SW, "_", corrname, ".jpg"),
         p1to4, width=9, height=7, units="in", dpi=800)
  
  return(p1to4)
}

varSCSW_line_plot_sequences <- function(S_vals, r, m_SW, reps_SW, m_SC, reps_SC,
                                   pre_SC, post_SC, corrtype, pereff, title=""){
  
  rhovals <- seq(0.01, 0.2, 0.005)
  r <- r
  relvars <- expand.grid(S=S_vals, rho=rhovals, r=r)
  relvars$varSC <- with(
    relvars,
    sapply(1:nrow(relvars), function(j){
      CRTVarSW(m_SC, SCdesmat(S[j], reps_SC, pre_SC, post_SC),
               rho[j], r[j], corrtype, pereff)
    }
    )
  )
  relvars$varSW <- with(
    relvars,
    sapply(1:nrow(relvars), function(j){
      CRTVarSW(m_SW, SWdesmat(S[j], reps_SW),
               rho[j], r[j], corrtype, pereff)
    }
    )
  )
  relvars <- relvars %>%
    mutate(relvarSCSW = varSC/varSW,
           Sfac=as.factor(S)) %>%
    select(c(rho, Sfac, relvarSCSW))
  
  p <- ggplot(data=relvars, aes(x=rho, y=relvarSCSW, colour=Sfac)) +
    geom_line(size=1.2) +
    geom_hline(yintercept=1, linetype="dashed") +
    #      expand_limits(y=ylimits) +
    xlab("Within-period ICC") +
    ylab("Relative variance") +
    labs(title=title, colour="Sequences") +
    theme_bw() +
    theme(plot.title=element_text(hjust=0.5, size=12),
          axis.title=element_text(size=10), axis.text=element_text(size=10),
          legend.key.width = unit(1.5, "cm"),
          legend.title=element_text(size=12), legend.text=element_text(size=12),
          legend.position="bottom")
  return(p)
}

varSCSW_multi_plot_sequences <- function(S_vals, r, reps_SW, reps_SC, pre_SC, post_SC, corrtype){
  # Compare variances of complete SW and staircase designs, for a range of
  # correlation parameters
  # Inputs:
  #  S_vals - numbers of unique treatment sequences, to form lines on plots, e.g. c(3, 5, 10, 20)
  #  r - cluster autocorrelation value for all scenarios
  #  reps_SW - number of times each sequence is repeated for SW design
  #  reps_SC - number of times each sequence is repeated for SC design
  #  pre_SC - number of pre-switch measurement periods for SC design
  #  post_SC - number of post-switch measurement periods for SC design
  #  corrtype - within-cluster correlation structure type
  #             (0=block-exchangeable, 1=exponential decay)
  # Output:
  #  Multiplot of relative variances (vartheta_SC/vartheta_SW), for m=10 and 100
  #  and categorical and linear period effects
  
  m1 <- 10
  m2 <- 100
  
  # m=10, categorical period effects
  p1 <- varSCSW_line_plot_sequences(
    S_vals, r,
    m1, reps_SW,
    m1, reps_SC, pre_SC, post_SC,
    corrtype, 'cat',
    bquote(paste("m = ", .(m1), ", ", "categorical period effects"))
  )
  # m=100, categorical period effects
  p2 <- varSCSW_line_plot_sequences(
    S_vals, r,
    m2, reps_SW,
    m2, reps_SC, pre_SC, post_SC,
    corrtype, 'cat',
    bquote(paste("m = ", .(m2), ", ", "categorical period effects"))
  )
  # m=10, linear period effects
  p3 <- varSCSW_line_plot_sequences(
    S_vals, r, 
    m1, reps_SW,
    m1, reps_SC, pre_SC, post_SC,
    corrtype, 'lin',
    bquote(paste("m = ", .(m1), ", ", "linear period effects"))
  )
  # m=100, linear period effects
  p4 <- varSCSW_line_plot_sequences(
    S_vals, r,
    m2, reps_SW,
    m2, reps_SC, pre_SC, post_SC,
    corrtype, 'lin',
    bquote(paste("m = ", .(m2), ", ", "linear period effects"))
  )
  mylegend <- g_legend(p1)
  title <- bquote(paste("Relative variance of treatment effect estimators, ",
                        Var(hat(theta))[paste("SC(S,", .(reps_SC), ",", .(pre_SC), ",", .(post_SC), ")")]/
                          Var(hat(theta))[paste("SW(S,", .(reps_SW), ")")]))
  p1to4 <- make_2x2_multiplot(p1, p2, p3, p4, mylegend, title=title)
  corrname <- ifelse(corrtype==0, "BE", "DTD")
  ggsave(paste0("plots/SC_S", reps_SC, pre_SC, post_SC, "_vs_SW_S",
                reps_SW, "_", corrname, ".jpg"),
         p1to4, width=9, height=7, units="in", dpi=800)
  return(p1to4)
}

VarSW_alt <- function(m, S, K, rho0, r){
  totalvar <- 1
  sig2CP <- rho0*totalvar
  sig2E <- totalvar - sig2CP
  sig2 <- sig2E/m
  a <- sig2 + sig2CP
  b <- r*sig2CP
  (12*(a - b)*(a + S*b))/(K*(S^2 - 1)*(2*a + S*b))
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

VarSCbasic_line_plot <- function(m, S, K, pereff){
  
  rhovals <- seq(0.01, .99, 0.01)
  rvals <- c(0.2, 0.5, 0.8, 1.0)
  vars <- expand.grid(rho=rhovals, r=rvals)
  
  if(pereff=="cat"){
    vars$varSC <- with(
      vars,
      sapply(1:nrow(vars), function(j){
        VarSCcat(m, S, K, rho[j], r[j]*rho[j], r[j])
      }
      )
    )
  }else if(pereff=="lin"){
    vars$varSC <- with(
      vars,
      sapply(1:nrow(vars), function(j){
        VarSClin(m, S, K, rho[j], r[j]*rho[j], r[j])
        }
      )
    )
  }

  vars <- vars %>%
    mutate(rfac=as.factor(r)) %>%
    select(c(rho, rfac, varSC))
  
  title <- bquote(paste("Variance of treatment effect estimator, ",
                        Var(hat(theta))[paste("SC(", .(S), ",", .(K), ",", "1,1),", .(pereff))]))
  p <- ggplot(data=vars, aes(x=rho, y=varSC, colour=rfac)) +
    geom_line(size=1.2) +
    xlab("Within-period ICC") +
    ylab("Variance") +
    labs(title=title, colour="Cluster autocorrelation") +
    theme_bw() +
    theme(plot.title=element_text(hjust=0.5, size=16),
          axis.title=element_text(size=14), axis.text=element_text(size=14),
          legend.key.width = unit(1.5, "cm"),
          legend.title=element_text(size=14), legend.text=element_text(size=14),
          legend.position="bottom")
  ggsave(paste0("plots/SC_", S, K, "11_", "m", m, "_", pereff, ".jpg"),
         p, width=9, height=7, units="in", dpi=800)
  return(p)
}

VarSCbasic_line_plot_m <- function(S, K, rho, pereff){
  
  mvals <- seq(10, 1000, 10)
  rvals <- c(0.2, 0.5, 0.8, 1.0)
  vars <- expand.grid(m=mvals, r=rvals)
  
  if(pereff=="cat"){
    vars <- vars %>%
      mutate(
        varSC = VarSCcat(m, S, K, rho, rho*r, r),
        S = S
      )
  }else if(pereff=="lin"){
    vars <- vars %>%
      mutate(
        varSC = VarSClin(m, S, K, rho, rho*r, r),
        S = S
      )
  }

  vars <- vars %>%
    mutate(rfac=as.factor(r)) %>%
    select(-r)
  
  title <- bquote(paste("Variance of treatment effect estimator, ",
                        Var(hat(theta))[paste("SC(", .(S), ",", .(K), ",", "1,1),", .(pereff))]))
  p <- ggplot(data=vars, aes(x=m, y=varSC, colour=rfac)) +
    geom_line(size=1.2) +
    xlab(expression(paste("Cluster-period size, ", m))) +
    ylab("Variance") +
    labs(title=title, colour=expression(paste("Cluster autocorrelation, ", r))) +
    theme_bw() +
    theme(plot.title=element_text(hjust=0.5, size=16),
          axis.title=element_text(size=16), axis.text=element_text(size=16),
          legend.key.width = unit(1.5, "cm"),
          legend.title=element_text(size=16), legend.text=element_text(size=16),
          legend.position="bottom")
  ggsave(paste0("plots/SC_", S, K, "11_rho_", rho, "_diffm_", pereff, ".jpg"),
         p, width=9, height=7, units="in", dpi=800)
  return(p)
}

pow <- function(vars, effsize, siglevel=0.05){
  z <- qnorm(siglevel/2)
  pow <- pnorm(z + sqrt(1/vars)*effsize)
  return(pow)
}
