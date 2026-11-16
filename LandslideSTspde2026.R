# Analysis for "Scalar-on-Function Regression with space-time SPDE Smoothing"
# Erin Bryce, Daniela Castro-Camilo, and David L Miller
# 2026

library(tidyr)
library(raster)
library(INLA)
library(ggplot2)
library(dplyr)
library(mgcv)
library(readr)
library(lubridate)

# functions for SPDE basis construction
source("smooth.construct.spdeST.smooth.spec.R")

## data

# raw data
# we have 23140 locations
# measured at 174 times each, 12 days apart

# date range we're interested in
start_date <- as.Date("2017-04-05")
end_date   <- as.Date("2018-11-25")

## get land deformation data 2017-04-05 - 2018-11-25
dat <- read_csv("data/dat.csv")
dat <- subset(dat, time>=start_date & time<=end_date)

# load static covariates
ogcovs <- read_csv("data/static features.csv")

# load precipitation data
# this has columns for dates and as many rows as we have
# location-observation times
precip <- read_csv("data/precip.csv")
# remove first/last columns
precip <- precip[, -c(1, ncol(precip))]
# dates available in precip data
precip_dates <- ymd(names(precip))
# convert to matrix
precip <- as.matrix(precip)

# number of days to use
ndays <- 50

# storage for precip data
precm <- matrix(NA, nrow(dat), ncol=ndays)

# to make this fast, avoid matching
base_date <- as.numeric(min(precip_dates))

# this takes a while, so create a progress bar
pb <- txtProgressBar(min=0, max=nrow(dat), style=3)

# create the distributed lag columns
for(i in 1:nrow(dat)){
  this_date <- dat[i,]$time
  this_site <- dat[i,]$su_id
  # index dates, most recent first: column 1 is the observation day,
  # column ndays is (ndays-1) days before it
  date_rng <- this_date:(this_date - days(ndays-1))

  precm[i, ] <- precip[this_site+1, date_rng-base_date+1]

  setTxtProgressBar(pb, i)
}
close(pb)


plot(c(1,ndays), range(precm), type="n")
for(i in sample(1:nrow(precm), 100)){
  lines(1:ndays, precm[i,])
}

# generate the lag matrix
# prect[, k] is the lag for precm[, k]: k=1 is the observation day (0 days
# before), k=ndays is ndays-1 days before. Lag is stored as days-before + 1 so
# it doubles as the column index (plot_dlnm_cumulative.R relies on this).
prect <- matrix(1:ndays, nrow(dat), ncol=ndays, byrow=TRUE)

## extract the relevant dates
#date_index <- as.Date(colnames(precip)) >= start_date &
#  as.Date(colnames(precip)) <= end_date
#
#precip <- precip[, date_index]
#
#chunks <- split(1:ncol(precip), ceiling(seq_along(1:ncol(precip))/12))
#chunks <- lapply(chunks, function(cols) precip[,cols])
#chunks <- lapply(chunks, function(df) {
#  colnames(df) <- paste0("prec", seq_len(ncol(df)))
#  df
#})
#
#precdf <- do.call(rbind, chunks)
## nrow = 1:23140 SU repeated 173 times and ncol = 12 days
#precm  <- as.matrix(precdf)
#
## create precpitation time matrix
#n_rows   <- 23140
#n_blocks <- 50
## create columns
## TODO: Dave check that this works?
##precl <- lapply(12:23, \(x){
##  rep((1:n_blocks - 1) * 1+(x/12), each = n_rows)
##})
### cbind them together into a data.frame
##prectime <- do.call(cbind, precl)
#prectime <- matrix(1:12, nrow=nrow(precm), ncol=12, byrow=TRUE)
#colnames(prectime) <- paste0("prec", 1:12)



# positive y values only
dat$sdeform <- dat$sdeform + 5

# fix data names
dat$lith   <- as.factor(dat$lithology)
dat$slope  <- dat$avg_slope
dat$faults <- dat$dis2faults
dat$rivers <- dat$dis2river
dat$planc  <- dat$avg_plan_c
dat$profc  <- dat$avg_prof_c



# distributed lag without refund
dat$precm <- precm
dat$prect  <- prect


# subset for quick fitting
#subind <- 1:231400
alldata <- dat#[subind, ]

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


b_dl <- bam(sdeform ~ te(precm, prect, k=20, bs="bs") +
                      s(slope, k=18, bs="bs") +
                      s(faults, k=18, bs="bs") +
                      s(rivers, k=18, bs="bs") +
                      lith +
                      planc +
                      profc +
                      s(sux, suy, timeID, bs = "spdeST",
                        xt = list(mesh = mesh.s, mesh.time = mesh.time)),
             data = alldata,
             family=Gamma(link = "log"),
             samfrac=0.2,
             # discretisation bins, one per smooth term in formula order.
             # 1D terms (and te() marginals) use mgcv's default of 1000. A 3D
             # term with m bins gets round(m^(1/3))+1 per dimension; the default
             # (25 per dimension) put the sites on a ~25x25 grid and merged
             # neighbouring times. 98^3 gives 99 per dimension: timeID 1..50
             # lands exactly on the grid (spacing 0.5) and sites move < 1 km.
             # Cost: ~420k unique rows x 8450 cols, about 28 GB dense
             discrete = c(1000, 1000, 1000, 1000, 98^3),
             nthreads = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4")),
             method = "fREML")
summary(b_dl)
# results identical to the one with refund, so just use this?

# save dl model
save(b_dl, file="dl_model.rda")


# save data
save(dat, file="data.rda")

