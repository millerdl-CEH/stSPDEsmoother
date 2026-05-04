# model comparison

library(ggplot2)
library(mgcv)
library(INLA)
source("smooth.construct.spdeST.smooth.spec.R")
library(gratia)

# load the fitted distributed lag model
load("dl_model.rda")
# SoFR model
load("sofr_model.rda")
load("data.rda")

# SoFR model seems to be better by AIC at least?
AIC(b_dl, b_sofr)

# not much difference between models, so interpretation is important?
plot(predict(b_dl), predict(b_sofr))



