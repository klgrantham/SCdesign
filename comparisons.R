# Efficiency comparisons for staircase designs
#
# Kelsey Grantham (kelsey.grantham@monash.edu)

source('variances.R')

## Block exchangeable within-cluster correlation ##

# SW design with 3 clusters, 4 periods vs.
# SC design with 3 clusters, 2 periods per sequence
# Note: This is the embedded staircase design
varSCSW_plot_general(10, 3, 1, 10, 3, 1, 1, 1, 0, 'cat')
varSCSW_plot_general(10, 3, 1, 10, 3, 1, 1, 1, 0, 'lin')

# SW design with 3 clusters, 4 periods vs.
# SC design with 6 clusters, 2 periods per sequence
# (6+1+1-1=7 periods total)
# Note: Designs have same total number of cluster-periods
varSCSW_plot_general(10, 3, 1, 10, 6, 1, 1, 1, 0, 'cat')
varSCSW_plot_general(10, 3, 1, 10, 6, 1, 1, 1, 0, 'lin')
# SC design with 6 clusters (3 sequences, each repeated twice), 2 periods each
varSCSW_plot_general(10, 3, 1, 10, 3, 1, 1, 2, 0, 'cat')
varSCSW_plot_general(10, 3, 1, 10, 3, 1, 1, 2, 0, 'lin')

# SW design with 3 clusters, 4 periods vs.
# SC design with 3 clusters, 4 periods per sequence
# (3+2+2-1=6 periods total)
# Note: Designs have same total number of cluster-periods
varSCSW_plot_general(10, 3, 1, 10, 3, 2, 2, 1, 0, 'cat')
varSCSW_plot_general(10, 3, 1, 10, 3, 2, 2, 1, 0, 'lin')

# SW design with 21 clusters (7 unique sequences), 8 periods vs.
# SC design with 21 clusters (7 unique sequences), 2 periods per sequence
# (7+1+1-1=8 periods total)
# Note: This is the embedded staircase design
varSCSW_plot_general(10, 7, 3, 10, 7, 1, 1, 3, 0, 'cat')
# Relative variance around 2 for most plausible correlation values in bottom-
# right quadrant (low ICC and high cluster autocorrelation)

# Relative variance is higher for designs with more cluster-periods, with the
# exception of some combinations with high ICC and high cluster autocorrelation

# 100 subjects per cluster-period compared to 10
varSCSW_plot_general(100, 7, 3, 100, 7, 1, 1, 3, 0, 'cat')
# Slightly higher relative variance for some configurations, lower for others

