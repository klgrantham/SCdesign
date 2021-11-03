# Efficiency comparisons for staircase designs
#
# Kelsey Grantham (kelsey.grantham@monash.edu)

source('variances.R')

library(ggplot2)
library(grid)
library(gridExtra)

# Extract legend
g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

# Multiplot framework
make_2x2_multiplot <- function(p1, p2, p3, p4, legend, title){
  p <- grid.arrange(arrangeGrob(p1 + theme(legend.position="none"),
                                p2 + theme(legend.position="none"),
                                p3 + theme(legend.position="none"),
                                p4 + theme(legend.position="none"),
                                ncol=2),
                    legend, nrow=2, heights=c(10,1),
                    top=textGrob(title,
                                 gp=gpar(fontsize=18)))
  return(p)
}

make_1x2_multiplot <- function(p1, p2, legend, title){
  p <- grid.arrange(arrangeGrob(p1 + theme(legend.position="none"),
                                p2 + theme(legend.position="none"),
                                ncol=2),
                    legend, nrow=2, heights=c(10,1),
                    top=textGrob(title,
                                 gp=gpar(fontsize=18)))
  return(p)
}

## Block exchangeable within-cluster correlation ##

# SW design with 3 clusters, 4 periods vs.
# SC design with 3 clusters, 2 periods per sequence
# Note: This is the embedded staircase design
varSCSW_grid_plot(10, 3, 1, 10, 3, 1, 1, 1, 0, 'cat')
varSCSW_grid_plot(10, 3, 1, 10, 3, 1, 1, 1, 0, 'lin')

# SW design with 3 clusters, 4 periods vs.
# SC design with 6 clusters, 2 periods per sequence
# (6+1+1-1=7 periods total)
# Note: Designs have same total number of cluster-periods
varSCSW_grid_plot(10, 3, 1, 10, 6, 1, 1, 1, 0, 'cat')
varSCSW_grid_plot(10, 3, 1, 10, 6, 1, 1, 1, 0, 'lin')
# SC design with 6 clusters (3 sequences, each repeated twice), 2 periods each
varSCSW_grid_plot(10, 3, 1, 10, 3, 2, 1, 1, 0, 'cat')
varSCSW_grid_plot(10, 3, 1, 10, 3, 2, 1, 1, 0, 'lin')

# SW design with 3 clusters, 4 periods vs.
# SC design with 3 clusters, 4 periods per sequence
# (3+2+2-1=6 periods total)
# Note: Designs have same total number of cluster-periods
varSCSW_grid_plot(10, 3, 1, 10, 3, 1, 2, 2, 0, 'cat')
varSCSW_grid_plot(10, 3, 1, 10, 3, 1, 2, 2, 0, 'lin')

# SW design with 21 clusters (7 unique sequences), 8 periods vs.
# SC design with 21 clusters (7 unique sequences), 2 periods per sequence
# (7+1+1-1=8 periods total)
# Note: This is the embedded staircase design
varSCSW_grid_plot(10, 7, 3, 10, 7, 3, 1, 1, 0, 'cat')
# Relative variance around 2 for most plausible correlation values in bottom-
# right quadrant (low ICC and high cluster autocorrelation)

# Relative variance is higher for designs with more cluster-periods, with the
# exception of some combinations with high ICC and high cluster autocorrelation

# 100 subjects per cluster-period compared to 10
varSCSW_grid_plot(100, 7, 3, 100, 7, 3, 1, 1, 0, 'cat')
# Slightly higher relative variance for some configurations, lower for others

varSCSW_line_plot(10, 3, 1, 10, 3, 1, 1, 1, 0, 'cat')

varSCSW_line_plot(10, 3, 1, 10, 3, 1, 1, 1, 0, 'lin')

varSCSW_line_plot(100, 3, 1, 100, 3, 1, 1, 1, 0, 'cat')
# For larger cluster-period size, relative variance ramps up more quickly as
# within-period ICC increases

# block-exchangeable
varSCSW_multi_plot(3, 1, 3, 1, 1, 1, 0)
# discrete-time decay
varSCSW_multi_plot(3, 1, 3, 1, 1, 1, 1)

varSCSW_line_plot(100, 7, 3, 100, 7, 3, 1, 1, 0, 'cat')
varSCSW_line_plot(100, 7, 1, 100, 7, 1, 1, 1, 0, 'cat')
# Relative variances remain the same if values of reps_SC and reps_SW change,
# as long as reps_SC'=w*reps_SC and reps_SW'=w*reps_SW, for some positive factor w

# m=100, block-exchangeable correlation, categorical period effects
varSCSW_line_plot(100, 3, 1, 100, 6, 1, 1, 1, 0, 'cat')
# For same number of cluster-periods, extended SC around twice as efficient as
# SW design for values of r<0.75

# block-exchangeable
varSCSW_multi_plot(3, 1, 6, 1, 1, 1, 0)
# discrete-time decay
varSCSW_multi_plot(3, 1, 6, 1, 1, 1, 1)

# m=100, discrete-time decaying correlation, categorical period effects
varSCSW_line_plot(100, 3, 1, 100, 6, 1, 1, 1, 'cat')

varSCSW_multi_plot(7, 5, 7, 5, 1, 1, 0)
varSCSW_multi_plot(7, 5, 7, 5, 1, 1, 1)
