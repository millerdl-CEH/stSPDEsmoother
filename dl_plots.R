# some code to make plots
library(ggplot2)
library(mgcv)
library(INLA)
source("smooth.construct.spdeST.smooth.spec.R")
library(gratia)

# load the fitted distributed lag model
load("dl_model.rda")
load("data.rda")


# plot some effects (parametric + smooths but not SPDE or distributed lag)
p_dl_eff <- draw(b_dl, select=c("s(slope)", "s(faults)", "s(rivers)"), parametric=TRUE) &
  theme_minimal()
p_dl_eff
ggsave(p_dl_eff, file="plots/dl_smooth_effects.pdf", width=10, height=8)


## plot SPDE

# first build the spatial bits
predsux   <- seq(min(alldata$sux), max(alldata$sux), length.out = 50)
predsuy   <- seq(min(alldata$suy), max(alldata$suy), length.out = 50)
predgrid  <- expand.grid(sux = predsux, suy = predsuy)

# which time indices do we want to show?
predtimes  <- c(1, 10, 25, 40, 50)
fitted_vals <- list()

# do this in bits to keep within vector size limits in R
for(pt in predtimes){
#  predgrid$timeID <- pt

  predgrid2 <- predgrid[c(rep(1, 50), 1:nrow(predgrid)),]
  # all group timeID values need to be in the prediction to get Xp to
  # be the right size
  predgrid2$timeID <- c(1:50, rep(pt, nrow(predgrid2)-50))

  Xp <- Predict.matrix.spdeST.smooth(b_dl$smooth[[5]], data = predgrid2)

  # remove those 50
  Xp <- Xp[-(1:50),]
# TODO: not clear why last column (all zeros) is here. It mismatches
#       the number of coefficients?
Xp <- Xp[,-ncol(Xp)]

  coefs        <- coef(b_dl)
  coefs_spde   <- coefs[grepl("s\\(sux,suy,timeID\\)", names(coefs))]
  fitted_vals[[pt]]  <- as.vector(Xp %*% coefs_spde)
}

# make a big prediction grid for plotting
predgrid <- predgrid[rep(1:nrow(predgrid), length(predtimes)), ]
predgrid$timeID <- predtimes[rep(1:length(predtimes), each=nrow(predgrid)/length((predtimes)))]
predgrid$fit <- do.call(c, fitted_vals)


p_dl_spde <- ggplot(predgrid, aes(x = sux, y = suy, fill = fit)) +
  geom_tile() +
  facet_wrap(~timeID) +
  scale_fill_viridis_c() +
  scale_x_continuous(expand=FALSE) +
  scale_y_continuous(expand=FALSE) +
  coord_equal() +
  theme_minimal() +
  labs(title = "SPDE smooth term")
p_dl_spde

ggsave(p_dl_spde, file="plots/dl_spde.pdf", width=10, height=6)

# distributed lag term

plot_grid <- expand.grid(precm = seq(min(alldata$precm),
                                     max(alldata$precm),
                                     length.out=100),
                         prect = seq(min(alldata$prect),
                                     max(alldata$prect),
                                     length.out=100),
                         slope=0,
                         faults=0,
                         rivers=0,
                         lith=alldata$lith[1],
                         planc=0,
                         profc=0,
                         sux=0,
                         suy=0,
                         timeID=1)

# make predictions per term
pr <- predict(b_dl, plot_grid, type="terms", terms="te(precm,prect)", discrete=FALSE)

# we just want the first column te(time, precip)
plot_grid$p <- pr[,1]

# make nice labels for the lag axis
#lag_breaks <- seq(0, 56, by=14)
#lag_labels <- format(ymd("1970-01-01") + lag_breaks, "%d %b")

# finally make the plot
p_dl_effect <- ggplot(plot_grid) +
  geom_tile(aes(x=precm, y=prect, fill=p)) +
  geom_contour(aes(x=precm, y=prect, z=p)) +
  metR::geom_text_contour(aes(x=precm, y=prect, z=p),
           rotate=FALSE, label.placer = metR::label_placer_flattest()) +
  scale_fill_gradient2(low = scales::muted("blue"),
                       high = scales::muted("red")) +
#  scale_y_continuous(breaks=lag_breaks, labels=lag_labels) +
  labs(x="Precipitation", y="Date", fill="Distributed\nlag effect") +
  theme_minimal() +
  coord_cartesian(expand=FALSE)
p_dl_effect

ggsave(p_dl_effect, file="plots/dl_effect.pdf", width=8, height=8)

## more DL plots

# build the Lp matrix, the design matrix for predictions
lp <- predict(b_dl, plot_grid, type="lpmatrix", terms="te(precm,prect)", discrete=FALSE)

# simulate from the posterior of the model
bs <- rmvn(1000, coef(b_dl), vcov(b_dl))
# zero the intercept
bs[, !grepl("te\\(precm,prect\\)", names(coef(b_dl)))] <- 0

# calculate predictions
pr <- lp%*%t(bs)

plot_grid$lower <- NA
plot_grid$upper <- NA
plot_grid$p <- NA

# for each precip
for(p in unique(plot_grid$precm)){
  # calculate the per-temperature, per-simulation effect
  precs <- apply(pr[plot_grid$precm==p, ], 2, sum)

  # now get summary statistics over the posterior samples
  qs <- quantile(precs, prob=c(0.025, 0.975))
  plot_grid[plot_grid$precm==p,]$lower <- qs[1]
  plot_grid[plot_grid$precm==p,]$upper <- qs[2]
  plot_grid[plot_grid$precm==p,]$p <- median(precs)
}

# make data for the rug plot
precrug <- data.frame(precm=c(alldata$precm))

# finally make the plot
p_prec_marginal <- ggplot(plot_grid) +
  geom_ribbon(aes(x=precm, ymin=lower, ymax=upper), alpha=0.6) +
  geom_line(aes(x=precm, y=p)) +
  geom_rug(aes(x=precm), data=precrug) +
  labs(x="Precipitation", y="Marginal effect") +
  theme_minimal()
p_prec_marginal

ggsave(p_prec_marginal, file="plots/dl_marginal.pdf", width=8, height=8)


