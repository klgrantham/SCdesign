# Variance and power calculations for a specific trial example
#
# Kelsey Grantham (kelsey.grantham@monash.edu)
#

source('variances.R')

# Trial design choices:
#   * Effect size of 0.5
#   * 50 subjects/cluster-period
#   * Exchangeable intracluster correlation (r=1)
#   * ICC of 0.1

### Basic staircase with S=4, K=1 ###
SC4111 <- SCdesmat(S=4, K=1, pre=1, post=1)

## Categorical period effects

# Using explicit variance expression for basic staircase designs
var_SC4111_cat <- VarSCcat(m=50, S=4, K=1, rho0=0.1, rhou=0.1, r=1)

# Using general variance expression
var_SC4111_cat <- CRTvartheta(m=50, SC4111, rho0=0.1, r=1, corrtype=0, pereff='cat')

pow_SC4111_cat <- pow(var_SC4111_cat, effsize=0.5, siglevel=0.05)

## Linear period effects

# Using explicit varince expression for basic staircase designs
var_SC4111_lin <- VarSClin(m=50, S=4, K=1, rho0=0.1, rhou=0.1, r=1)

# Using general variance expression
var_SC4111_lin <- CRTvartheta(m=50, SC4111, rho0=0.1, r=1, corrtype=0, pereff='lin')

pow_SC4111_lin <- pow(var_SC4111_lin, effsize=0.5, siglevel=0.05)

### Imbalanced staircase with S=4, K=1, R0=1, R1=3 ###
SC4113 <- SCdesmat(S=4, K=1, pre=1, post=3)

## Categorical period effects

# Using general variance expression
var_SC4113_cat <- CRTvartheta(m=50, SC4113, rho0=0.1, r=1, corrtype=0, pereff='cat')
pow_SC4113_cat <- pow(var_SC4113_cat, effsize=0.5, siglevel=0.05)

## Linear period effects

# Using general variance expression
var_SC4113_lin <- CRTvartheta(m=50, SC4113, rho0=0.1, r=1, corrtype=0, pereff='lin')
pow_SC4113_lin <- pow(var_SC4113_lin, effsize=0.5, siglevel=0.05)
