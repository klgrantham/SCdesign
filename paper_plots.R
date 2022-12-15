# Variance plots for staircase design
#
# Kelsey Grantham (kelsey.grantham@monash.edu)

source('variances.R')

Varvals_SCbasic_rho <- function(m, S, K, pereff){
  
  rhovals <- seq(0.01, .99, 0.01)
  rvals <- c(0.2, 0.5, 0.8, 1.0)
  vars <- expand.grid(rho=rhovals, r=rvals)
  
  if(pereff=="cat"){
    vars <- vars %>%
      mutate(
        varSC = VarSCcat(m, S, K, rho, r*rho, r),
        m = m,
        S = S
      )
  }else if(pereff=="lin"){
    vars <- vars %>%
      mutate(
        varSC = VarSClin(m, S, K, rho, r*rho, r),
        m = m,
        S = S
      )
  }
  return(vars)  
}

VarSCbasic_multi_line_plot <- function(pereff){
  # Compare variances of basic staircase designs for different configurations
  # for a range of correlation parameters
  # Inputs:
  #  pereff - time period effect type
  #           ('cat'=categorical period effects, 'lin'=linear period effects)
  # Output:
  #  Multiplot of variances (vartheta_SC), for m=10 and 100 (rows),
  #    for S=3 and S=10 (columns) for the specified time effect
  
  S1 <- 3
  S2 <- 10
  m1 <- 10
  m2 <- 100
  
  vars_m1_S1 <- Varvals_SCbasic_rho(m1, S1, 1, pereff)
  vars_m1_S2 <- Varvals_SCbasic_rho(m1, S2, 1, pereff)
  vars_m2_S1 <- Varvals_SCbasic_rho(m2, S1, 1, pereff)
  vars_m2_S2 <- Varvals_SCbasic_rho(m2, S2, 1, pereff)

  allvars <- bind_rows(
    vars_m1_S1,
    vars_m1_S2,
    vars_m2_S1,
    vars_m2_S2
  )
  
  allvars <- allvars %>%
    mutate(rfac = as.factor(r)) %>%
    select(-r)
  
  m.labs <- c("m = 10", "m = 100")
  names(m.labs) <- c(m1, m2)
  S.labs <- c("S = 3", "S = 10")
  names(S.labs) <- c(S1, S2)
  
  title <- bquote(paste("Variance of treatment effect estimator, ",
                        Var(hat(theta))[paste("SC(S,1,1,1),", .(pereff))]))
  p <- ggplot(allvars, aes(x=rho, y=varSC, colour=rfac, linetype=rfac)) +
    geom_line(size=1.2) +
    facet_grid(
      m ~ S,
      labeller = labeller(m = m.labs, S = S.labs)
    ) +
    xlab(expression(paste("Within-period ICC, ", rho[0]))) +
    ylab("Variance") +
    labs(title=title, colour=expression(paste("Cluster autocorrelation, ", r)),
         linetype=expression(paste("Cluster autocorrelation, ", r))) +
    theme_bw() +
    theme(plot.title=element_text(hjust=0.5, size=16),
          axis.title=element_text(size=14), axis.text=element_text(size=14),
          strip.background = element_rect(
            color="white", fill="white", linetype="solid"
          ),
          strip.text.x = element_text(size = 14),
          strip.text.y = element_text(size=14),
          legend.key.width = unit(1.5, "cm"),
          legend.title=element_text(size=14), legend.text=element_text(size=14),
          legend.position="bottom")
  ggsave(paste0("plots/multiplot_SCbasic_S_", S1, "vs", S2, "_m_", m1, "vs", m2, "_", pereff, ".jpg"),
         p, width=9, height=5, units="in", dpi=300)
  return(p)
}


Varvals_SCbasic_m <- function(rho, S, K, pereff){
  
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
  return(vars)  
}

VarSCbasic_line_plot_m <- function(rho, S, K, pereff, ylims=NA, compare=FALSE, refvars=NA){
  
  vars <- Varvals_SCbasic_m(rho, S, K, pereff)
  
  vars <- vars %>%
    mutate(rfac=as.factor(r)) %>%
    select(-r)
  
  title <- bquote(paste("Variance of treatment effect estimator, ",
                        Var(hat(theta))[paste("SC(", .(S), ",", .(K), ",", "1,1),", .(pereff))]))
  p <- ggplot(data=vars, aes(x=m, y=varSC, colour=rfac)) +
    geom_line(size=1.2) +
    ylim(ylims) +
    xlab(expression(paste("Cluster-period size, ", m))) +
    ylab("Variance") +
    labs(title=title, colour=expression(paste("Cluster autocorrelation, ", r))) +
    theme_bw() +
    theme(plot.title=element_text(hjust=0.5, size=18),
          axis.title=element_text(size=16), axis.text=element_text(size=16),
          legend.key.width = unit(1.5, "cm"),
          legend.title=element_text(size=16), legend.text=element_text(size=16),
          legend.position="bottom")
  
  if(compare==TRUE){
    refvars <- refvars %>%
      mutate(rfac=as.factor(r)) %>%
      select(-r)
    
    p <- p +
      geom_line(data=refvars, size=1.2, alpha=0.2, aes(x=m, y=varSC, colour=rfac))
  }
  
  ggsave(paste0("plots/SC_", S, K, "11_rho_", rho, "_diffm_", pereff, ".jpg"),
         p, width=9, height=7, units="in", dpi=800)
  
  return(p)
}

VarSCbasic_multi_line_plot_m <- function(rhos, S, K, pereff, ylims=c(0,1)){
  
  rho1 <- rhos[1]
  rho2 <- rhos[2]

  vars_rho1 <- Varvals_SCbasic_m(rho1, S, K, pereff)
  vars_rho1$rho <- rho1
  vars_rho2 <- Varvals_SCbasic_m(rho2, S, K, pereff)
  vars_rho2$rho <- rho2
  
  vars <- bind_rows(
    vars_rho1,
    vars_rho2
  )
  
  vars <- vars %>%
    mutate(rfac=as.factor(r),
           rhofac=as.factor(rho)) %>%
    select(-r, -rho)

  title <- bquote(paste("Variance of treatment effect estimator, ",
                        Var(hat(theta))[paste("SC(", .(S), ",", .(K), ",", "1,1),", .(pereff))]))
  p <- ggplot(data=vars, aes(x=m, y=varSC, colour=rfac, linetype=rfac)) +
    geom_line(size=1.2) +
    facet_grid(
      . ~ rhofac,
      labeller = label_bquote(cols = rho[0]==.(as.character(rhofac)))
    ) +
    ylim(ylims) +
    xlab(expression(paste("Cluster-period size, ", m))) +
    ylab("Variance") +
    labs(title=title, colour=expression(paste("Cluster autocorrelation, ", r)),
         linetype=expression(paste("Cluster autocorrelation, ", r))) +
    theme_bw() +
    theme(plot.title=element_text(hjust=0.5, size=16),
          axis.title=element_text(size=14), axis.text=element_text(size=14),
          strip.background = element_rect(
            color="white", fill="white", linetype="solid"
          ),
          strip.text.x = element_text(size = 14),
          strip.text.y = element_text(size=14),
          legend.key.width = unit(1.5, "cm"),
          legend.title=element_text(size=14), legend.text=element_text(size=14),
          legend.position="bottom")

  ggsave(paste0("plots/SC_", S, K, "11_rhos_", rho1, "_", rho2, "_diffm_", pereff, ".jpg"),
         p, width=9, height=5, units="in", dpi=300)
  
  return(p)
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
  
  rhovals <- seq(0.01, 0.25, 0.01)
  rvals <- seq(0.2, 0.95, 0.05)
  vars <- expand.grid(rho=rhovals, r=rvals)
  
  vars$varSCmat <- with(
    vars,
    sapply(1:nrow(vars), function(j){
      CRTVarSW(m_SC, SCdesmat(S_SC, reps_SC, pre_SC, post_SC),
               rho[j], r[j], corrtype, pereff)
      }
    )
  )
  vars$varSW <- with(
    vars,
    sapply(1:nrow(vars), function(j){
      CRTVarSW(m_SW, SWdesmat(S_SW, reps_SW),
               rho[j], r[j], corrtype, pereff)
      }
    )
  )
  
  if(pereff=="cat"){
    vars <- vars %>%
      mutate(
        varSC = VarSCcat(m_SC, S_SC, reps_SC, rho, r*rho, r),
      )
  }else if(pereff=="lin"){
    vars <- vars %>%
      mutate(
        varSC = VarSClin(m_SC, S_SC, reps_SC, rho, r*rho, r),
      )
  }
  
  vars <- vars %>%
    mutate(
      relvarSCSW = varSC/varSW,
      relvarSCSWmat = varSCmat/varSW,
      releffSCSW = varSW/varSC
    )
  return(vars)
}


releffSCSW_grid_multiplot_corr <- function(S_SW, reps_SW, S_SC, reps_SC,
                                           pre_SC, post_SC, pereff,
                                           fixedscale=FALSE, limits=c(1,2),
                                           breaks=seq(1,2,0.5)){
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
  relvars_m1_BE <- gridvals_small(m1, S_SW, reps_SW, m1, S_SC, reps_SC,
                                  pre_SC, post_SC, 0, pereff)
  relvars_m1_BE$m <- "m = 10"
  relvars_m1_BE$corrname <- "Block-exchangeable"
  relvars_m2_BE <- gridvals_small(m2, S_SW, reps_SW, m2, S_SC, reps_SC,
                                  pre_SC, post_SC, 0, pereff)
  relvars_m2_BE$m <- "m = 100"
  relvars_m2_BE$corrname <- "Block-exchangeable"
  relvars_m1_DTD <- gridvals_small(m1, S_SW, reps_SW, m1, S_SC, reps_SC,
                                   pre_SC, post_SC, 1, pereff)
  relvars_m1_DTD$m <- "m = 10"
  relvars_m1_DTD$corrname <- "Discrete-time decay"
  relvars_m2_DTD <- gridvals_small(m2, S_SW, reps_SW, m2, S_SC, reps_SC,
                                   pre_SC, post_SC, 1, pereff)
  relvars_m2_DTD$m <- "m = 100"
  relvars_m2_DTD$corrname <- "Discrete-time decay"
  relvars <- bind_rows(
    relvars_m1_BE,
    relvars_m2_BE,
    relvars_m1_DTD,
    relvars_m2_DTD
  )
  
  if(fixedscale==TRUE){
    fillopt <- scale_fill_viridis_c(name="Relative efficiency", direction=-1,
                                    limits=limits, breaks=breaks)
  }else{
    fillopt <- scale_fill_viridis_c(name="Relative efficiency", direction=-1)
  }
  
  p <- ggplot(relvars, aes(x=r, y=rho)) +
    geom_tile(aes(fill=releffSCSW)) +
    fillopt +
    facet_grid(
      m ~ corrname
    ) +
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
    ggtitle(bquote(paste("Relative efficiency, ",
                         Var(hat(theta))[paste("SW(", .(S_SW), ",", .(reps_SW), ")")]/
                           Var(hat(theta))[paste("SC(", .(S_SC), ",", .(reps_SC), ",", .(pre_SC), ",", .(post_SC), ")")])))
  ggsave(paste0("plots/releff_SC_", S_SC, reps_SC, pre_SC, post_SC, "_vs_SW_",
                S_SW, reps_SW, "_", pereff, ".jpg"),
         p, width=9, height=5, units="in", dpi=800)
  return(p)
}

releffSCSW_grid_multiplot_diffm <- function(S_SW, reps_SW, S_SC, reps_SC,
                                            pre_SC, post_SC, pereff,
                                            fixedscale=FALSE, limits=c(1,2),
                                            breaks=seq(1,2,0.5)){
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
  
  m1_SW <- 10
  m1_SC <- ((S_SW+1)/2)*m1_SW
  m2_SW <- 100
  m2_SC <- ((S_SW+1)/2)*m2_SW
  relvars_m1_BE <- gridvals_small(m1_SW, S_SW, reps_SW, m1_SC, S_SC, reps_SC,
                                  pre_SC, post_SC, 0, pereff)
  relvars_m1_BE$m <- "m1"
  relvars_m1_BE$corrname <- "Block-exchangeable"
  relvars_m2_BE <- gridvals_small(m2_SW, S_SW, reps_SW, m2_SC, S_SC, reps_SC,
                                  pre_SC, post_SC, 0, pereff)
  relvars_m2_BE$m <- "m2"
  relvars_m2_BE$corrname <- "Block-exchangeable"
  relvars_m1_DTD <- gridvals_small(m1_SW, S_SW, reps_SW, m1_SC, S_SC, reps_SC,
                                   pre_SC, post_SC, 1, pereff)
  relvars_m1_DTD$m <- "m1"
  relvars_m1_DTD$corrname <- "Discrete-time decay"
  relvars_m2_DTD <- gridvals_small(m2_SW, S_SW, reps_SW, m2_SC, S_SC, reps_SC,
                                   pre_SC, post_SC, 1, pereff)
  relvars_m2_DTD$m <- "m2"
  relvars_m2_DTD$corrname <- "Discrete-time decay"
  relvars <- bind_rows(
    relvars_m1_BE,
    relvars_m2_BE,
    relvars_m1_DTD,
    relvars_m2_DTD
  )
  
  if(fixedscale==TRUE){
    fillopt <- scale_fill_viridis_c(name="Relative efficiency", direction=-1,
                                    limits=limits, breaks=breaks)
  }else{
    fillopt <- scale_fill_viridis_c(name="Relative efficiency", direction=-1)
  }
  
  m.labs <- c(paste("mSC = ", m1_SC, "\nmSW = 10"), paste("mSC = ", m2_SC, "\nmSW = 100"))
#  m.labs <- c("mSC=(S+1)mSW/2,\nmSW=10", "mSC=(S+1)mSW/2,\nmSW=100")
  names(m.labs) <- c("m1", "m2")
  
  p <- ggplot(relvars, aes(x=r, y=rho)) + 
    geom_tile(aes(fill=releffSCSW)) +
    facet_grid(
      m ~ corrname,
      labeller = labeller(m = m.labs)
    ) +
    fillopt +
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
          strip.text.x = element_text(size = 10),
          strip.text.y = element_text(size=10)) +
    coord_fixed() + xlab(expression(paste("Cluster autocorrelation, ", r))) + ylab(expression(paste("Within-period ICC, ", rho))) +
    ggtitle(bquote(paste("Relative efficiency, ",
                         Var(hat(theta))[paste("SW(", .(S_SW), ",", .(reps_SW), ")")]/
                           Var(hat(theta))[paste("SC(", .(S_SC), ",", .(reps_SC), ",", .(pre_SC), ",", .(post_SC), ")")])))
  ggsave(paste0("plots/releff_diffm_SC_", S_SC, reps_SC, pre_SC, post_SC, "_vs_SW_",
                S_SW, reps_SW, "_", pereff, ".jpg"),
         p, width=9, height=5, units="in", dpi=800)
  return(p)
}

# Variance of treatment effect estimator, varying within-period ICC
VarSCbasic_multi_line_plot('cat')
VarSCbasic_multi_line_plot('lin')

# Variance of treatment effect estimator, varying cluster-period size
VarSCbasic_line_plot_m(0.05, 10, 1, 'cat', ylims=c(0, 0.035))
VarSCbasic_line_plot_m(0.05, 10, 1, 'lin', ylims=c(0, 0.035))

VarSCbasic_multi_line_plot_m(c(0.05, 0.2), 10, 1, 'cat', ylims=c(0, 0.06))
VarSCbasic_multi_line_plot_m(c(0.05, 0.2), 10, 1, 'lin', ylims=c(0, 0.06))

## Relative efficiency

# Embedded staircase vs stepped wedge
releffSCSW_grid_multiplot_corr(3, 1, 3, 1, 1, 1, 'cat', fixedscale=TRUE, limits=c(0,1), breaks=seq(0,1,0.2))
releffSCSW_grid_multiplot_corr(3, 1, 3, 1, 1, 1, 'lin', fixedscale=TRUE, limits=c(0,1), breaks=seq(0,1,0.2))

releffSCSW_grid_multiplot_corr(10, 1, 10, 1, 1, 1, 'cat', fixedscale=TRUE, limits=c(0,1), breaks=seq(0,1,0.2))
releffSCSW_grid_multiplot_corr(10, 1, 10, 1, 1, 1, 'lin', fixedscale=TRUE, limits=c(0,1), breaks=seq(0,1,0.2))

# Staircase with larger cluster-period size vs stepped wedge
releffSCSW_grid_multiplot_diffm(3, 1, 3, 1, 1, 1, 'cat')
releffSCSW_grid_multiplot_diffm(3, 1, 3, 1, 1, 1, 'cat', fixedscale=TRUE, limits=c(0,2.5), breaks=seq(0,2.5,0.5))
releffSCSW_grid_multiplot_diffm(3, 1, 3, 1, 1, 1, 'lin', fixedscale=TRUE, limits=c(0,2.5), breaks=seq(0,2.5,0.5))

releffSCSW_grid_multiplot_diffm(10, 1, 10, 1, 1, 1, 'cat', fixedscale=TRUE, limits=c(0,2.5), breaks=seq(0,2.5,0.5))
releffSCSW_grid_multiplot_diffm(10, 1, 10, 1, 1, 1, 'lin', fixedscale=TRUE, limits=c(0,2.5), breaks=seq(0,2.5,0.5))

# Extended staircase vs stepped wedge
releffSCSW_grid_multiplot_corr(3, 1, 3, 2, 1, 1, 'cat', fixedscale=TRUE, limits=c(0,2.0), breaks=seq(0,2.0,0.5))
releffSCSW_grid_multiplot_corr(3, 1, 3, 2, 1, 1, 'lin', fixedscale=TRUE, limits=c(0,2.0), breaks=seq(0,2.0,0.5))


## Trial examples

# PROMPT trial
CRTVarSW(20, SWdesmat(5, 1), 0.03, 1, 0, 'cat')/VarSCcat(20, 5, 1, 0.03, 0.03, 1)
CRTVarSW(20, SWdesmat(5, 1), 0.03, 1, 0, 'cat')/VarSCcat(33, 5, 1, 0.03, 0.03, 1)
