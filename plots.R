# some code to make plots
library(ggplot2)
library(mgcv)
library(INLA)
source("smooth.construct.spdeST.smooth.spec.R")
library(gratia)

# load the fitted distributed lag model
load("dl_model.rda")
b_dl <- model
load("data.rda")
dat <- subset


# plot some effects (parametric + smooths but not SPDE or distributed lag)
draw(b_dl, select=c("s(slope)", "s(faults)", "s(rivers)"), parametric=TRUE) &
  theme_minimal()

# plot SPDE


# distributed lag term

plot_grid <- expand.grid(precm.tmat = seq(min(dat$precm.tmat),
                                          max(dat$precm.tmat),
                                          length.out=100),
                         precm.omat = seq(min(dat$precm.omat),
                                          max(dat$precm.omat),
                                          length.out=100),
                         slope=0,
                         faults=0,
                         rivers=0,
                         lith=dat$lith[1],
                         planc=0,
                         profc=0,
                         timeID=1)

plot_grid$sux <- seq(min(dat$sux), max(dat$sux), length.out=nrow(plot_grid))
plot_grid$suy <- seq(min(dat$suy), max(dat$suy), length.out=nrow(plot_grid))

#plot_grid <- dat

# make predictions per term
pr <- predict(b_dl, plot_grid, type="terms", terms="te(precm.tmat,precm.omat)", discrete=FALSE)

# we just want the first column s(lag, temp)
plot_grid$p <- pr[,1]

# make nice labels for the lag axis
#lag_breaks <- seq(0, 56, by=14)
#lag_labels <- format(ymd("1970-01-01") + lag_breaks, "%d %b")

# finally make the plot
p_dl_effect <- ggplot(plot_grid) +
  geom_tile(aes(x=precm.tmat, y=precm.omat, fill=p)) +
  geom_contour(aes(x=precm.tmat, y=precm.omat, z=p)) +
  metR::geom_text_contour(aes(x=precm.tmat, y=precm.omat, z=p),
           rotate=FALSE, label.placer = metR::label_placer_flattest()) +
  scale_fill_gradient2(low = scales::muted("blue"),
                       high = scales::muted("red")) +
#  scale_y_continuous(breaks=lag_breaks, labels=lag_labels) +
  labs(x="Precipitation", y="Date", fill="Distributed\nlag effect") +
  theme_minimal() +
  coord_cartesian(expand=FALSE)
p_dl_effect

