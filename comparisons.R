# Efficiency comparisons for staircase designs
#
# Kelsey Grantham (kelsey.grantham@monash.edu)

source('variances.R')

# SW design with 3 clusters, 4 periods vs.
# SC design with 3 clusters, 2 periods per sequence
# Note: This is the embedded staircase design
varSCSW_plot_general(10, 3, 1, 10, 3, 1, 1, 1, 0, 'cat')

# SW design with 3 clusters, 4 periods vs.
# SC design with 6 clusters, 2 periods per sequence
# (6+1+1-1=7 periods total)
# Note: Designs have same total number of cluster-periods
varSCSW_plot_general(10, 3, 1, 10, 6, 1, 1, 1, 0, 'cat')
# SC design with 6 clusters (3 sequences, each repeated twice), 2 periods each
varSCSW_plot_general(10, 3, 1, 10, 3, 1, 1, 2, 0, 'cat')

# SW design with 3 clusters, 4 periods vs.
# SC design with 3 clusters, 4 periods per sequence
# (3+2+2-1=6 periods total)
# Note: Designs have same total number of cluster-periods
varSCSW_plot_general(10, 3, 1, 10, 3, 2, 2, 1, 0, 'cat')
varSCSW_plot_general(10, 3, 1, 10, 3, 2, 2, 1, 0, 'lin')
