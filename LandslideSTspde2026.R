# Analysis for "Scalar-on-Function Regression with space-time SPDE Smoothing"
# Erin Bryce, Daniela Castro-Camilo, and David L Miller
# 2026

library(tidyr)
library(raster)
library(INLA)
library(spdep)
library(ggplot2)
library(inlabru)
library(dplyr)
library(mgcv)
library(mgcViz)
library(readr)

# functions for SPDE basis construction
source("smooth.construct.spdeST.smooth.spec.R")

#### DATA ####
ogcovs <- read_csv("data/static features.csv")

## get data 2017-04-16 - 2018-11-16
dat <- read_csv("data/dat.csv")

dat$time <- as.POSIXct(dat$time, format = "%Y-%m-%d")

precip <- read_csv("data/precip.csv")
precip <- precip[,-c(1, 2193)]
colnames(precip) <- as.Date(sub("^X", "", names(precip)), format = "%Y%m%d")
start_date <- as.Date("2017-04-05")
end_date   <- as.Date("2018-11-25")
date_index <- as.Date(colnames(precip)) >= start_date &
  as.Date(colnames(precip)) <= end_date

precip <- precip[, date_index]

chunks <- split(1:ncol(precip), ceiling(seq_along(1:ncol(precip))/12))
chunks <- lapply(chunks, function(cols) precip[,cols])
chunks <- lapply(chunks, function(df) {
  colnames(df) <- paste0("prec", seq_len(ncol(df)))
  df
})

precdf <- do.call(rbind, chunks)
# nrow = 1:23140 SU repeated 173 times and ncol = 12 days
precm  <- as.matrix(precdf)

# create precpitation time matrix
n_rows   <- 23140
n_blocks <- 50
# create columns
precl <- lapply(12:23, \(x){
  rep((1:n_blocks - 1) * 1+(x/12), each = n_rows)
})
# cbind them together into a data.frame
prectime <- do.call(cbind, precl)
colnames(prectime) <- paste0("prec", 1:12)


dat$time <- as.Date(dat$time)
data     <- dat %>%
  filter(time >= start_date & time <= end_date)

# positive y values only
data$sdeform <- data$sdeform + 5
alldata      <- cbind(data, precdf)

#### Matern penalty in space-time ####

mesh.s <- inla.mesh.2d(loc = cbind(alldata$sux, alldata$suy),
                       max.edge = .2*c(1, 2),
                       cutoff = .1)

# Setting up the Q_time matrix
T         <- length(unique(alldata$timeID))
delta_t   <- 1
a         <- 0.5
I         <- Diagonal(T)
Q_time    <- (1 / delta_t^2) * I + (a^2) * I
mesh.time <- list(qt = Q_time)

# update ar1 matrix:
T       <- length(unique(alldata$timeID))
delta_t <- 1
a       <- 0.5
phi     <- exp(-a * delta_t)

main_diag <- c(1, rep(1 + phi^2, T-2), 1)
off_diag  <- rep(-phi, T-1)
Q_time <- bandSparse(
  T, T,
  k = c(-1, 0, 1),
  diagonals = list(off_diag, main_diag, off_diag)
)
Q_time <- Q_time / (1 - phi^2)
mesh.time <- list(qt = Q_time)

alldata$lith   <- as.factor(alldata[, "lithology"])

## distributed lag model with af
#library(refund)
#paf <- af(precm, argvals = prectime)
#alldata$precm.tmat <- paf$data$precm.tmat
#alldata$precm.omat <- paf$data$precm.omat
#
#b_dlr <- bam(sdeform ~ te(precm.tmat, precm.omat) +
#                       s(slope, k=15) +
#                       s(faults, k=15) +
#                       s(rivers) +
#                       lith +
#                       planc +
#                       profc +
#                       s(sux, suy, timeID, bs = "spdeST",
#                         xt = list(mesh = mesh.s, mesh.time = mesh.time)),
#             data = alldata,
#             family=Gamma(link = "log"),
#             discrete = TRUE, nthreads = 4,
#             method = "fREML")
#summary(b_dlr)
#
## save dl model
#save(b_dlr, file="dlr_model.rda")

# distributed lag without refund
alldata$precm <- precm
alldata$prect  <- prectime

b_dl <- bam(sdeform ~ te(precm, prect) +
                      s(slope, k=15) +
                      s(faults, k=15) +
                      s(rivers) +
                      lith +
                      planc +
                      profc +
                      s(sux, suy, timeID, bs = "spdeST",
                        xt = list(mesh = mesh.s, mesh.time = mesh.time)),
             data = alldata,
             family=Gamma(link = "log"),
             discrete = TRUE, nthreads = 4,
             method = "fREML")
summary(b_dl)
# results identical to the one with refund, so just use this?

# save dl model
save(b_dl, file="dl_model.rda")


# scalar-on-function regression model

b_sofr <- bam(sdeform ~ s(by=precm, prect) +
                      s(slope, k=15) +
                      s(faults, k=15) +
                      s(rivers) +
                      lith +
                      planc +
                      profc +
                      s(sux, suy, timeID, bs = "spdeST",
                        xt = list(mesh = mesh.s, mesh.time = mesh.time)),
             data = alldata,
             family=Gamma(link = "log"),
             discrete = TRUE, nthreads = 4,
             method = "fREML")
summary(b_sofr)

# save dl model
save(b_sofr, file="sofr_model.rda")

# SoFR model seems to be better by AIC at least?
AIC(b_dlr, b_dl, b_sofr)

# not much difference between models, so interpretation is important?
plot(predict(b_dl), predict(b_sofr))

# save data
save(alldata, file="data.rda")

