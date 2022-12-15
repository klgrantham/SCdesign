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
                        var(hat(theta))[paste("SC(S,1,1,1),", .(pereff))]))
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

## Keep ##
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
                        var(hat(theta))[paste("SC(", .(S), ",", .(K), ",", "1,1),", .(pereff))]))
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
  return(p)
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

# Figure 2: Variance of treatment effect estimator, varying within-period ICC,
#           categorical period effects
p2 <- VarSCbasic_multi_line_plot('cat')
ggsave("plots/multiplot_SCbasic_S_3vs10_m_10vs100_cat.jpg",
       p2, width=9, height=5, units="in", dpi=300)
ggsave(paste0("plots/figure2.eps"), p2, width=9, height=5, units="in", dpi=800)

# Figure 3: Variance of treatment effect estimator, varying cluster-period size,
#           categorical period effects
p3 <- VarSCbasic_multi_line_plot_m(c(0.05, 0.2), 10, 1, 'cat', ylims=c(0, 0.06))
ggsave(paste0("plots/SCbasic_S_10_K_1_rhos_5_20_diffm_cat.jpg"),
       p3, width=9, height=5, units="in", dpi=300)
ggsave(paste0("plots/figure3.eps"), p3, width=9, height=5, units="in", dpi=800)

# Figure 4: Variance of treatment effect estimator, varying within-period ICC,
#           linear period effects
p4 <- VarSCbasic_multi_line_plot('lin')
ggsave(paste0("plots/multiplot_SCbasic_S_3vs10_m_10vs100_lin.jpg"),
       p4, width=9, height=5, units="in", dpi=300)
ggsave(paste0("plots/figure4.eps"), p4, width=9, height=5, units="in", dpi=800)

# Figure 5: Variance of treatment effect estimator, varying cluster-period size,
#           linear period effects
p5 <- VarSCbasic_multi_line_plot_m(c(0.05, 0.2), 10, 1, 'lin', ylims=c(0, 0.06))
ggsave(paste0("plots/SCbasic_S_10_K_1_rhos_5_20_diffm_lin.jpg"),
       p5, width=9, height=5, units="in", dpi=300)
ggsave(paste0("plots/figure5.eps"), p5, width=9, height=5, units="in", dpi=800)
