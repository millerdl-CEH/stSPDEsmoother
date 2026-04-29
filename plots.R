# some code to make plots
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
dat <- subset


# plot some effects (parametric + smooths but not SPDE or distributed lag)
draw(b_dl, select=c("s(slope)", "s(faults)", "s(rivers)"), parametric=TRUE) &
  theme_minimal()

# plot SPDE
predsux   <- seq(min(subset$sux), max(subset$sux), length.out = 50)
predsuy   <- seq(min(subset$suy), max(subset$suy), length.out = 50)
predtime  <- seq(min(subset$timeID), max(subset$timeID), by = 1)
predgrid  <- expand.grid(sux = predsux, suy = predsuy, timeID = predtime)  # ideally whole {x,y}
Xp <- Predict.matrix.spdeST.smooth(b_dl$smooth[[5]], data = predgrid)

coefs        <- coef(b_dl)
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

# distributed lag term

plot_grid <- expand.grid(precm = seq(min(dat$precm),
                                     max(dat$precm),
                                     length.out=100),
                         prect = seq(min(dat$prect),
                                     max(dat$prect),
                                     length.out=100),
                         slope=0,
                         faults=0,
                         rivers=0,
                         lith=dat$lith[1],
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

ggsave(p_dl_effect, file="dl_effect.pdf", width=8, height=8)

## more DL plots

# build the Lp matrix, the design matrix for predictions
lp <- predict(b_dl, plot_grid, type="lpmatrix", terms="te(precm,prect)", discrete=FALSE)

# simulate from the posterior of the model
# TODO: switch to mh sampling
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
precrug <- data.frame(precm=c(dat$precm))

# finally make the plot
p_prec_marginal <- ggplot(plot_grid) +
  geom_ribbon(aes(x=precm, ymin=lower, ymax=upper), alpha=0.6) +
  geom_line(aes(x=precm, y=p)) +
  geom_rug(aes(x=precm), data=precrug) +
  labs(x="Precipitation", y="Marginal effect") +
  theme_minimal()
p_prec_marginal

ggsave(p_prec_marginal, file="dl_marginal.pdf", width=8, height=8)




### stuff to move


qqplot(subset$sdeform, fitted(b_dl), xlim = c(0, 60), ylim = c(0, 60))
abline(0, 1)


## plots lithology and curvature
#theme_plot <- function() {
#  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
#        panel.grid.minor = element_blank(),
#        legend.position = "none",
#        plot.title = element_blank(),
#        axis.text.y = element_text(size=12),
#        axis.text.x = element_text(size=12),
#        axis.title.x = element_blank(),
#        plot.margin = unit(c(1,3,1,1), "lines"),
#        axis.ticks.length.x = unit(.3, "cm"))
#}
#gridExtra::grid.arrange(
#  dplyr::filter(fe, group == "curvature") |>
#    ggplot(aes(x = name, y = rc))+
#    geom_errorbar(aes(ymin = rc-2*se, ymax = rc+2*se),
#                  linewidth=0.75, width = 0.5) +
#    geom_point(aes(x = name, y = rc), size = 2, color = 'red') +
#    theme_plot(),
#
#  dplyr::filter(fe, group == "lithology") |>
#    ggplot(aes(x = name, y = rc))+
#    geom_errorbar(aes(ymin = rc-2*se, ymax = rc+2*se),
#                  linewidth=0.75, width = 0.5) +
#    geom_point(aes(x = name, y = rc), size = 2, color = 'red') +
#    theme_plot(),
#
#  ncol = 2)


