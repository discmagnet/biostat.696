

#################
#################    Case-control studies: test on ratio of intensity function
#################    and Cuzick-Edwards test
#################

library(MASS)
library(maptools)
library(splancs)
library(rgdal)
library(spatialkernel)
library(sp)
library(spatstat)
library(spdep)
library(fields)
library(SpatialEpi)
library(smacpod)

data(grave)
## Plotting the grave data: location of 143 graves, of which some of are affected
## and some are not affected
plot(grave)

## Computing the ratio of the intensity function for the affected and the non affected ## graves
renv = logrr(grave)
plot(renv)

## To perform a hypothesis test, we need to proceed via simulations and carry out
## a Monte Carlo test.
## Here we perform the Monte Carlo test, running 50 permutations of the label
renv = logrr(grave, nsim = 50)
logrr.test(renv)
[1] "The p-value for the global test is 0.196078431372549"

## Cuzick-Edwards test
qnn.test(grave,q=c(3,5,7,9,11,13,15))
$qsum
q  Tq pvalue
1  3  32  0.004
2  5  45  0.022
3  7  58  0.040
4  9  73  0.040
5 11  91  0.014
6 13 109  0.008
7 15 122  0.010

$consum
contrast Tcontrast pvaluecon
1    T5 - T3        13     0.422
2    T7 - T3        26     0.376
3    T9 - T3        41     0.256
4   T11 - T3        59     0.114
5   T13 - T3        77     0.050
6   T15 - T3        90     0.064
7    T7 - T5        13     0.454
8    T9 - T5        28     0.252
9   T11 - T5        46     0.086
10  T13 - T5        64     0.032
11  T15 - T5        77     0.038
12   T9 - T7        15     0.266
13  T11 - T7        33     0.056
14  T13 - T7        51     0.028
15  T15 - T7        64     0.034
16  T11 - T9        18     0.054
17  T13 - T9        36     0.022
18  T15 - T9        49     0.046
19 T13 - T11        18     0.102
20 T15 - T11        31     0.138
21 T15 - T13        13     0.418


