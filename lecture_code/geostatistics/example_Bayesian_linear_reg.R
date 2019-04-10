
## Example of a Bayesian linear regression model using MCMCpack

library(MCMCpack)

## Dataset on standardized measures of fertility and socio economic indicators for 47 French-speaking provinces in Switzerland
data(swiss)

## The MCMCregress command runs a linear regression model with Gaussian error; inverse Gamma prior on the variance of the residuals, and multivariate Gaussian prior on the regression coefficients.
## mcmc specifies the number of MCMC iterations
## burnin specifies the number of MCMC iterations to discard for burn-in
## thin indicates the thinning
## beta.start provides the initial values for the regression coefficients. Default value is NA
## b0 is the prior mean for the regression coefficients. b0 should be a vector of the same length as the numbers of coefficients. If b0 is a scalar, then all the regression coefficients have the same prior mean.
## B0 is the inverse of the prior covariance matrix in the multivariate Gaussian prior for the beta. B0 should be a square matrix with number of rows and columns equal to the number of regression coefficients. If B0 is a scalar, then the inverse of the covariance matrix (and the covariance matrix) is taken to be the identity times the scalar entered for B0. Default value for B0 is 0, assuming independence among the regression coefficients a priori and an infinite variance (e.g. improper uniform prior on each regression coefficient)
## c0: is twice the shape parameter for the inverse Gamma prior on the variance of the residuals (sigma^2)
## d0" is twice the scale parameter for the inverse Gamma prior on the variance of the residuals (sigma^2)
## sigma.mu: alternatively, the prior on the the variance of the residuals can be specified so that the Inverse Gamma prior have shape and scale parameters so that sigma.mu is the prior mean for the variance of the residuals.
## sigma.var: alternatively, the prior on the the variance of the residuals can be specified so that the Inverse Gamma prior have shape and scale parameters so that sigma.var is the prior variance for the variance of the residuals.
## verbose: if 0 no message is provided on the display indicating at what point is the MCMC algorithm, if it's a number it provides display messages every number of iterations.

rm(list=ls())
bayes.lm <- MCMCregress(Fertility ~ Agriculture + Examination +                 Education + Catholic + Infant.Mortality,data=swiss,b0=c(0.5,0,0,0,0,0),sigma.mu=0.1,sigma.var=0.1,mcmc=10000,burnin=5000,verbose=500)

summary(bayes.lm)
plot(bayes.lm)

## Running Geweke diagnostic
geweke.diag(bayes.lm)
Fraction in 1st window = 0.1
Fraction in 2nd window = 0.5

(Intercept)      Agriculture      Examination        Education         Catholic Infant.Mortality
-0.4659          -1.2827           0.7403          -1.8690           2.2138           0.9247
sigma2
-0.1398


## Looking at autocorrelation function plot
acfplot(bayes.lm)

## Or doing them by hand
par(mfrow=c(2,2))
acf(bayes.lm[,2],main="ACF plot for beta - Agriculture")
acf(bayes.lm[,3],main="ACF plot for beta - Examination")
acf(bayes.lm[,4],main="ACF plot for beta - Education")
acf(bayes.lm[,7],main="ACF plot for sigma^2")

