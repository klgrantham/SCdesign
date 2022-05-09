# Variance and relative variance plots for staircase design
# Presentation plots
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

VarSCbasic_line_plot <- function(m, S, K, pereff, ylims=NA, compare=FALSE, refvars=NA){

  vars <- Varvals_SCbasic_rho(m, S, K, pereff)
  
  vars <- vars %>%
    mutate(rfac=as.factor(r)) %>%
    select(-r)
  
  title <- bquote(paste("Variance of treatment effect estimator, ",
                        Var(hat(theta))[paste("SC(", .(S), ",", .(K), ",", "1,1),", .(pereff))]))
  p <- ggplot(data=vars, aes(x=rho, y=varSC, colour=rfac)) +
    geom_line(size=1.2) +
    ylim(ylims) +
    xlab(expression(paste("Within-period ICC, ", rho))) +
    ylab("Variance") +
    labs(title=title, colour=expression(paste("Cluster autocorrelation, ", r))) +
    theme_bw() +
    theme(plot.title=element_text(hjust=0.5, size=16),
          axis.title=element_text(size=16), axis.text=element_text(size=16),
          legend.key.width = unit(1.5, "cm"),
          legend.title=element_text(size=16), legend.text=element_text(size=16),
          legend.position="bottom")
  
  if(compare==TRUE){
    refvars <- refvars %>%
      mutate(rfac=as.factor(r)) %>%
      select(-r)
    
    p <- p +
      geom_line(data=refvars, size=1.2, alpha=0.2, aes(x=rho, y=varSC, colour=rfac))
  }

  ggsave(paste0("plots/presentation/SC_", S, K, "11_", "m", m, "_", pereff, ".jpg"),
         p, width=9, height=7, units="in", dpi=800)
  
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
    theme(plot.title=element_text(hjust=0.5, size=16),
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
  
  ggsave(paste0("plots/presentation/SC_", S, K, "11_rho_", rho, "_diffm_", pereff, ".jpg"),
         p, width=9, height=7, units="in", dpi=800)
  
  return(p)
}

varSCSW_grid_plot <- function(m_SW, S_SW, reps_SW, m_SC, S_SC, reps_SC, pre_SC,
                              post_SC, corrtype, pereff, fixedscale=FALSE, limits=c(1,2), breaks=seq(1,2,0.5)){
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
  
  relvars <- gridvals_small(m_SW, S_SW, reps_SW, m_SC, S_SC, reps_SC,
                            pre_SC, post_SC, corrtype, pereff)
  
  if(fixedscale==TRUE){
    fillopt <- scale_fill_viridis_c(name="Relative variance", direction=-1,
                                    limits=limits, breaks=breaks)
  }else{
    fillopt <- scale_fill_viridis_c(name="Relative variance", direction=-1)
  }

  p <- ggplot(relvars, aes(x=rseq, y=rhoseq)) + 
    geom_tile(aes(fill=value)) +
    fillopt +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    theme(aspect.ratio=3/8,
          legend.key.width = unit(1.5, "cm"),
          legend.title=element_text(size=16), legend.text=element_text(size=16),
          legend.position="bottom",
          plot.title=element_text(hjust=0.5, size=16),
          axis.title=element_text(size=16), axis.text=element_text(size=16)) +
    coord_fixed() + xlab(expression(paste("Cluster autocorrelation, ", r))) + ylab(expression(paste("Within-period ICC, ", rho))) +
    ggtitle(bquote(paste("Relative variance of treatment effect estimators, ",
                         Var(hat(theta))[paste("SC(", .(S_SC), ",", .(reps_SC), ",", .(pre_SC), ",", .(post_SC), ")")]/
                           Var(hat(theta))[paste("SW(", .(S_SW), ",", .(reps_SW), ")")])))
  corrname <- ifelse(corrtype==0, "BE", "DTD")
  ggsave(paste0("plots/presentation/multiplot_SC_", S_SC, reps_SC, pre_SC, post_SC, "_vs_SW_",
                S_SW, reps_SW, "_", corrname, "_", pereff, ".jpg"),
         p, width=9, height=5, units="in", dpi=800)
  return(p)
}

releffSCSW_grid_plot <- function(m_SW, S_SW, reps_SW, m_SC, S_SC, reps_SC, pre_SC,
                              post_SC, corrtype, pereff, fixedscale=FALSE, limits=c(1,2), breaks=seq(1,2,0.5)){
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
  
  relvars <- gridvals_small(m_SW, S_SW, reps_SW, m_SC, S_SC, reps_SC,
                            pre_SC, post_SC, corrtype, pereff)
  
  relvars <- relvars %>%
    mutate(releff = 1/value)
  
  if(fixedscale==TRUE){
    fillopt <- scale_fill_viridis_c(name="Relative efficiency", direction=-1,
                                    limits=limits, breaks=breaks)
  }else{
    fillopt <- scale_fill_viridis_c(name="Relative efficiency", direction=-1)
  }
  
  p <- ggplot(relvars, aes(x=rseq, y=rhoseq)) + 
    geom_tile(aes(fill=releff)) +
    fillopt +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    theme(aspect.ratio=3/8,
          legend.key.width = unit(1.5, "cm"),
          legend.title=element_text(size=16), legend.text=element_text(size=16),
          legend.position="bottom",
          plot.title=element_text(hjust=0.5, size=16),
          axis.title=element_text(size=16), axis.text=element_text(size=16)) +
    coord_fixed() + xlab(expression(paste("Cluster autocorrelation, ", r))) + ylab(expression(paste("Within-period ICC, ", rho))) +
    ggtitle(bquote(paste("Relative efficiency, ",
                         (1/Var(hat(theta))[paste("SC(", .(S_SC), ",", .(reps_SC), ",", .(pre_SC), ",", .(post_SC), ")")])/
                           (1/Var(hat(theta))[paste("SW(", .(S_SW), ",", .(reps_SW), ")")]))))
  corrname <- ifelse(corrtype==0, "BE", "DTD")
  ggsave(paste0("plots/presentation/releff_SC_", S_SC, reps_SC, pre_SC, post_SC, "_m", m_SC, "_vs_SW_",
                S_SW, reps_SW, "_m", m_SW, "_", corrname, "_", pereff, ".jpg"),
         p, width=9, height=5, units="in", dpi=800)
  return(p)
}

# Variance plot: S=3, m=10
VarSCbasic_line_plot(10, 3, 1, 'lin', ylims=c(0, 0.85))

# Variance plot: S=10, m=10 (with S=3, m=10 overlaid)
VarSCbasic_line_plot(10, 10, 1, 'lin', ylims=c(0, 0.85), compare=TRUE, refvars=Varvals_SCbasic(10, 3, 1, 'lin'))

# Variance plot: S=3, m=100 (with S=3, m=10 overlaid)
VarSCbasic_line_plot(100, 3, 1, 'lin', ylims=c(0, 0.85), compare=TRUE, refvars=Varvals_SCbasic(10, 3, 1, 'lin'))
VarSCbasic_line_plot(100, 10, 1, 'lin', ylims=c(0, 0.2), compare=TRUE, refvars=Varvals_SCbasic(10, 10, 1, 'lin'))

# Variance plot: S=10, rho=0.05, varying m
VarSCbasic_line_plot_m(0.05, 10, 1, 'lin', ylims=c(0, 0.03))
VarSCbasic_line_plot_m(0.2, 10, 1, 'lin', ylims=c(0, 0.05), compare=TRUE, refvars=Varvals_SCbasic_m(0.05, 10, 1, 'lin'))

# Variance plot: S=10, m=10, categorical versus linear
VarSCbasic_line_plot(10, 10, 1, 'cat', ylims=c(0, 0.2), compare=TRUE, refvars=Varvals_SCbasic(10, 10, 1, 'lin'))

### Relative variance results: for discrete-time decay model, categorical period effects?

## Embedded staircase vs complete stepped wedge

# Relative variance plot: Contour plot for S=3, m=10
varSCSW_grid_plot(10, 3, 1, 10, 3, 1, 1, 1, 1, 'cat')
varSCSW_grid_plot(10, 3, 1, 10, 3, 1, 1, 1, 1, 'cat', fixedscale=TRUE, limits=c(1,6.5), breaks=seq(1,6.5,1.0))
releffSCSW_grid_plot(10, 3, 1, 10, 3, 1, 1, 1, 1, 'cat', fixedscale=TRUE, limits=c(0,1), breaks=seq(0,1,0.2))

# Relative variance plot: Contour plot for S=10, m=10 (present with similar legend shading?)
varSCSW_grid_plot(10, 10, 1, 10, 10, 1, 1, 1, 1, 'cat')
varSCSW_grid_plot(10, 10, 1, 10, 10, 1, 1, 1, 1, 'cat', fixedscale=TRUE, limits=c(1,6.5), breaks=seq(1,6.5,1.0))
releffSCSW_grid_plot(10, 10, 1, 10, 10, 1, 1, 1, 1, 'cat', fixedscale=TRUE, limits=c(0,1), breaks=seq(0,1,0.2))

# Relative variance plot: Contour plot for S=3, m=100
varSCSW_grid_plot(100, 3, 1, 100, 3, 1, 1, 1, 1, 'cat')
varSCSW_grid_plot(100, 3, 1, 100, 3, 1, 1, 1, 1, 'cat', fixedscale=TRUE, limits=c(1,6.5), breaks=seq(1,6.5,1.0))
releffSCSW_grid_plot(100, 3, 1, 100, 3, 1, 1, 1, 1, 'cat', fixedscale=TRUE, limits=c(0,1), breaks=seq(0,1,0.2))

# Relative variance plot: Contour plot for S=10, m=100
varSCSW_grid_plot(100, 10, 1, 100, 10, 1, 1, 1, 1, 'cat')
varSCSW_grid_plot(100, 10, 1, 100, 10, 1, 1, 1, 1, 'cat', fixedscale=TRUE, limits=c(1,6.5), breaks=seq(1,6.5,1.0))
releffSCSW_grid_plot(100, 10, 1, 100, 10, 1, 1, 1, 1, 'cat', fixedscale=TRUE, limits=c(0,1), breaks=seq(0,1,0.2))

## Relative efficiency: Staircase with larger cluster-period size vs stepped wedge

# Contour plot for S_SC=3, m_SC=20, S_SW=3, m_SW=10
releffSCSW_grid_plot(10, 3, 1, 20, 3, 1, 1, 1, 1, 'cat')
releffSCSW_grid_plot(10, 3, 1, 20, 3, 1, 1, 1, 1, 'cat', fixedscale=TRUE, limits=c(0.5,1.6), breaks=seq(0.5,1.5,0.25))
releffSCSW_grid_plot(10, 3, 1, 20, 3, 1, 1, 1, 1, 'cat', fixedscale=TRUE, limits=c(0,1.6), breaks=seq(0,1.5,0.5))

# Contour plot for S_SC=3, m_SC=200, S_SW=3, m_SW=100
releffSCSW_grid_plot(100, 3, 1, 200, 3, 1, 1, 1, 1, 'cat')
releffSCSW_grid_plot(100, 3, 1, 200, 3, 1, 1, 1, 1, 'cat', fixedscale=TRUE, limits=c(0,2.2), breaks=seq(0,2.0,0.5))
releffSCSW_grid_plot(100, 3, 1, 200, 3, 1, 1, 1, 1, 'cat', fixedscale=TRUE, limits=c(0,1.6), breaks=seq(0,1.5,0.5))

# Contour plot for S_SC=10, m_SC=55, S_SW=10, m_SW=10
releffSCSW_grid_plot(10, 10, 1, 55, 10, 1, 1, 1, 1, 'cat')
releffSCSW_grid_plot(10, 10, 1, 55, 10, 1, 1, 1, 1, 'cat', fixedscale=TRUE, limits=c(0,2.5), breaks=seq(0,2.5,0.5))

# Contour plot for S_SC=10, m_SC=550, S_SW=10, m_SW=100
releffSCSW_grid_plot(100, 10, 1, 550, 10, 1, 1, 1, 1, 'cat')
releffSCSW_grid_plot(100, 10, 1, 550, 10, 1, 1, 1, 1, 'cat', fixedscale=TRUE, limits=c(0,2.5), breaks=seq(0,2.5,0.5))

## Relative efficiency: Staircase with more clusters vs stepped wedge

# Contour plot for S_SC=6, m_SC=10, S_SW=3, m_SW=10
releffSCSW_grid_plot(10, 3, 1, 10, 6, 1, 1, 1, 1, 'cat')

# Contour plot for S_SC=3, K_SC=2, m_SC=10, S_SW=3, m_SW=10
releffSCSW_grid_plot(10, 3, 1, 10, 3, 2, 1, 1, 1, 'cat')
