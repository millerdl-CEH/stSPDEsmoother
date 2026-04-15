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
library(refund)
library(mgcv)
library(mgcViz)
library(readr)

# functions for SPDE basis construction
source("smooth.construct.spdeST.smooth.spec.R")

#### DATA ####
ogcovs <- read_csv("data/static features.csv")

## get data 2017-04-16 - 2018-11-16
dat <- read_csv("data/dat.csv")

precip <- read_csv("data/precip.csv")
precip <- precip[,-c(1, 2193)]
colnames(precip) <- as.Date(sub("^X", "", names(precip)), format = "%Y%m%d")
start_date <- as.Date("2017-04-16")
end_date   <- as.Date("2018-11-24")
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
n_blocks <- 49
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

subind <- 1:231400
subset <- alldata[subind, ]
mesh.s <- inla.mesh.2d(loc = cbind(subset$sux, subset$suy),
                       max.edge = .2*c(1, 2),
                       cutoff = .1)

# Setting up the Q_time matrix
T         <- length(unique(subset$timeID))
delta_t   <- 1
a         <- 0.5
I         <- Diagonal(T)
Q_time    <- (1 / delta_t^2) * I + (a^2) * I
mesh.time <- list(qt = Q_time)

# update ar1 matrix:
T       <- length(unique(subset$timeID))
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


subset$slope  <- as.vector(alldata[subind, "avg_slope"])
subset$faults <- as.vector(alldata[subind, "dis2faults"])
subset$rivers <- as.vector(alldata[subind, "dis2river"])
subset$lith   <- as.vector(as.factor(alldata[subind, "lithology"]))
subset$planc  <- as.vector(alldata[subind, "avg_plan_c"])
subset$profc  <- as.vector(alldata[subind, "avg_prof_c"])

precm <- precm[subind,]
prect  <- prectime[subind,]

paf <- af(precm, argvals = prect)
subset$precm.tmat <- paf$data$precm.tmat
subset$precm.omat <- paf$data$precm.omat
subset$L.precm    <- paf$data$L.precm


model <- bam(sdeform ~ te(precm.tmat, precm.omat) +
                       s(slope, k=15) +
                       s(faults) +
                       s(rivers) +
                       as.factor(lith) +
                       planc +
                       profc +
                       s(sux, suy, timeID, bs = "spdeST",
                         xt = list(mesh = mesh.s, mesh.time = mesh.time)),
             data = subset,
             family=Gamma(link = "log"),
             discrete = TRUE, nthreads = 4,
             method = "fREML")
summary(model)

predsux   <- seq(min(subset$sux), max(subset$sux), length.out = 50)
predsuy   <- seq(min(subset$suy), max(subset$suy), length.out = 50)
predtime  <- seq(min(subset$timeID), max(subset$timeID), by = 1)
predgrid  <- expand.grid(sux = predsux, suy = predsuy, timeID = predtime)  # ideally whole {x,y}
Xp <- Predict.matrix.spde.smooth(model$smooth[[5]], data = predgrid)

coefs        <- coef(model)
coefs_spde   <- tail(coefs, ncol(Xp))
fitted_vals  <- as.vector(Xp %*% coefs_spde)
predgrid$fit <- fitted_vals

ggplot(predgrid, aes(x = sux, y = suy, fill = fit)) +
  geom_tile() +
  facet_wrap(~timeID,
             labeller = labeller(timeID = setNames(paste0("Time point ",
               1:12), 1:12))) +
  scale_fill_viridis_c() +
  theme_classic() +
  labs(title = "SPDE smooth term")

qqplot(subset$sdeform, fitted(model), xlim = c(0, 60), ylim = c(0, 60))
abline(0, 1)

b <- getViz(model)
plot(b, select=1, too.far=0, main = "Precipitation")
#plot(model, select=1, scheme=2, too.far=0, main = "Precipitation") + l_dens(type = "cond") + l_fitLine() + l_ciLine()
plot(model, select=2, ylab="", xlab="Slope")
plot(model, select=3, ylab="", xlab="Distance to Faults")
plot(model, select=4, ylab="", xlab="Distance to Rivers")

fe <- dplyr::tibble(name = c("L2", "L4", "L5", "L8", "L10", "L11",
                             "planc", "profc"),
                    group = c(rep("lithology", 6), "curvature", "curvature"),
                    rc = c(summary(model)$p.coeff[4:7],
                           summary(model)$p.coeff[2]
                           summary(model)$p.coeff[3]
                           summary(model)$p.coeff[8:9]),
                    se = c(summary(model)$p.coeff[4:7]
                           summary(model)$p.coeff[2]
                           summary(model)$p.coeff[3]
                           summary(model)$p.coeff[8:9])) |>
  dplyr::mutate_at(c('name', 'group'), as.factor)

# plots lithology and curvature
theme_plot <- function() {
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        plot.title = element_blank(),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12),
        axis.title.x = element_blank(),
        plot.margin = unit(c(1,3,1,1), "lines"),
        axis.ticks.length.x = unit(.3, "cm"))
}
gridExtra::grid.arrange(
  dplyr::filter(fe, group == "curvature") |>
    ggplot(aes(x = name, y = rc))+
    geom_errorbar(aes(ymin = rc-2*se, ymax = rc+2*se),
                  linewidth=0.75, width = 0.5) +
    geom_point(aes(x = name, y = rc), size = 2, color = 'red') +
    theme_plot(),

  dplyr::filter(fe, group == "lithology") |>
    ggplot(aes(x = name, y = rc))+
    geom_errorbar(aes(ymin = rc-2*se, ymax = rc+2*se),
                  linewidth=0.75, width = 0.5) +
    geom_point(aes(x = name, y = rc), size = 2, color = 'red') +
    theme_plot(),

  ncol = 2)


