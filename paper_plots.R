# Variance plots for staircase design
#
# Kelsey Grantham (kelsey.grantham@monash.edu)

source('variances.R')

# Variance of treatment effect estimator, varying different parameters
VarSCbasic_line_plot(10, 3, 1, 'cat')
VarSCbasic_line_plot(10, 3, 1, 'lin')

VarSCbasic_line_plot(10, 10, 1, 'cat')
VarSCbasic_line_plot(10, 10, 1, 'lin')

VarSCbasic_multi_line_plot('cat')
VarSCbasic_multi_line_plot('lin')