
### This code is to fit a logistic regression model on individual level
### outcomes for which we don't have the exact geographical coordinates, however we
### do know in which areal unit the individuals are located.
###
### The data we have is of the type:
###  - we have N areal units
###  - in each areal unit, we have m_1, m_2,.., m_N subjects
###  - we have two covariates X1 and X2 on the individual
###  - we have a covariate Z on the areal unit where the individual lives
###  - we want to account for spatial correlation at the areal unit level
###
### We envision the following model:
###
###  log(odds for subject i in areal unit j)= beta_0 + beta_1 X1_i + beta_2 X2_i + beta_3 Z_j +
###                                           eta_j + epsilon_ij
###
### We model the spatial correlation in the spatial random effects eta_j using a CAR model
### For this we know the adjacency matrix for the areal units

## This code is to create the W matrix. The W matrix has 1 when the distance between two areal
## units in the square have a distance of 1 (meaning they share a border)
x.easting <- 1:10
x.northing <- 1:10
Grid <- expand.grid(x.easting, x.northing)
K <- nrow(Grid)

#### set up distance and neighbourhood (W, based on sharing a common border) matrices
distance <- as.matrix(dist(Grid))
W <-array(0, c(K,K))
W[distance==1] <-1


## Reading the data in
binary.data.ex <- read.table("Binary_data_multilevel_ex.txt",sep=" ",header=T)
names(binary.data.ex)

Y <- binary.data.ex$Y
X1 <- binary.data.ex$subj_covar1
X2 <- binary.data.ex$subj_covar2
Z <- binary.data.ex$area_covar
area <- binary.data.ex$area_subj
## Since the data is binary
trials <- rep(1,length(Y))

formula <- Y ~ X1 + X2 + Z
multi.model <- S.CARmultilevel(formula=formula, family="binomial", ind.area=area,
trials=trials, W=W, burnin=20000, fix.rho=T, rho=1,n.sample=100000)

names(multi.model)
[1] "summary.results"     "samples"             "fitted.values"       "residuals"           "modelfit"
[6] "accept"              "localised.structure" "formula"             "model"               "X"


multi.model$summary.results
             Median    2.5%   97.5% n.sample % accept n.effective Geweke.diag
(Intercept)  0.9615  0.5125  1.4021    80000     46.2      8006.0         0.2
X1           0.4757  0.3517  0.6027    80000     46.2     18291.0        -0.3
X2           0.0157 -0.1026  0.1348    80000     46.2     26927.2         0.2
Z           -0.8538 -1.0275 -0.6817    80000     46.2      8453.7        -0.1
tau2         0.1977  0.0612  0.4972    80000    100.0      1043.9        -1.0
rho          1.0000  1.0000  1.0000       NA       NA          NA          NA

## To make plots, you can proceed as done in the other code.






