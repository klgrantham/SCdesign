# Variance and relative variance plots for staircase design
# Presentation plots
#
# Kelsey Grantham (kelsey.grantham@monash.edu)

source('variances.R')

Varvals_SCbasic <- function(m, S, K, pereff){

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

  vars <- Varvals_SCbasic(m, S, K, pereff)
  
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
    labs(title=title, colour="Cluster autocorrelation") +
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

# Variance plot: S=3, m=10
VarSCbasic_line_plot(10, 3, 1, 'lin', ylims=c(0, 0.85))

# Variance plot: S=10, m=10 (with S=3, m=10 overlaid)
VarSCbasic_line_plot(10, 10, 1, 'lin', ylims=c(0, 0.85), compare=TRUE, refvars=Varvals_SCbasic(10, 3, 1, 'lin'))

# Variance plot: S=3, m=100 (with S=3, m=10 overlaid)
VarSCbasic_line_plot(100, 3, 1, 'lin', ylims=c(0, 0.85), compare=TRUE, refvars=Varvals_SCbasic(10, 3, 1, 'lin'))
VarSCbasic_line_plot(100, 10, 1, 'lin', ylims=c(0, 0.2), compare=TRUE, refvars=Varvals_SCbasic(10, 10, 1, 'lin'))

## Relative variance results: for discrete-time decay model, categorical period effects?

# Relative variance plot: Contour plot for S=3, m=10

# Relative variance plot: Contour plot for S=10, m=10 (present with similar legend shading?)

# Relative variance plot: Contour plot for S=3, m=100

# Relative variance plot: Contour plot for S=10, m=100
