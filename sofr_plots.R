# some code to make plots
library(ggplot2)
library(mgcv)
library(INLA)
source("smooth.construct.spdeST.smooth.spec.R")
library(gratia)
library(dplyr)
library(tidyr)
library(patchwork)

# SoFR model
load("sofr_model.rda")
load("data.rda")


# plot some effects (parametric + smooths but not SPDE or SoFR)
p_sofr_eff <- draw(b_sofr, select=c("s(slope)", "s(faults)", "s(rivers)"), parametric=TRUE) &
  theme_minimal()
p_sofr_eff
ggsave(p_sofr_eff, file="plots/sofr_smooth_effects.pdf", width=10, height=8)


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

  Xp <- Predict.matrix.spdeST.smooth(b_sofr$smooth[[5]], data = predgrid2)

  # remove those 50
  Xp <- Xp[-(1:50),]
# TODO: not clear why last column (all zeros) is here. It mismatches
#       the number of coefficients?
Xp <- Xp[,-ncol(Xp)]

  coefs        <- coef(b_sofr)
  coefs_spde   <- coefs[grepl("s\\(sux,suy,timeID\\)", names(coefs))]
  fitted_vals[[pt]]  <- as.vector(Xp %*% coefs_spde)
}

# make a big prediction grid for plotting
predgrid <- predgrid[rep(1:nrow(predgrid), length(predtimes)), ]
predgrid$timeID <- predtimes[rep(1:length(predtimes), each=nrow(predgrid)/length((predtimes)))]
predgrid$fit <- do.call(c, fitted_vals)


p_sofr_spde <- ggplot(predgrid, aes(x = sux, y = suy, fill = fit)) +
  geom_tile() +
  facet_wrap(~timeID) +
  scale_fill_viridis_c() +
  scale_x_continuous(expand=FALSE) +
  scale_y_continuous(expand=FALSE) +
  coord_equal() +
  theme_minimal() +
  labs(title = "SPDE smooth term")
p_sofr_spde

ggsave(p_sofr_spde, file="plots/sofr_spde.pdf", width=10, height=6)

# SoFR term

# setup the data
pp <- data.frame(prect=seq(1, 50, length.out=200),
                 precm=1,
                 slope=0,
                 faults=0,
                 rivers=0,
                 lith=alldata$lith[1],
                 planc=0,
                 profc=0,
                 sux=0,
                 suy=0,
                 timeID=1)

# make the prediction, per term and return standard errors
pr <- predict(b_sofr, type="terms", newdata=pp, se=TRUE, discrete=FALSE, terms="s(prect):precm")

# store the smooth term and calculate the confidence interval
pp$pr <- pr$fit
pp$ub <- pr$fit + 1.96*pr$se.fit
pp$lb <- pr$fit - 1.96*pr$se.fit

# construct the plot
precm_term <- ggplot(pp) +
  geom_ribbon(aes(x=prect, ymin=lb, ymax=ub), fill="grey70") +
  geom_line(aes(x=prect, y=pr)) +
  # zero line for comparison
  geom_hline(yintercept=0, colour="grey50", lty=2) +
#  labs(x="Day",
       # extract the EDF from the summary output
#       y=paste0("s(Day, ", round(summary(m_c)$edf[2], 2) ,")")) +
  # give informative labels
  scale_x_continuous(expand=expansion(0, 0)) +
  theme_minimal()
precm_term

rind <- sample(1:nrow(alldata), 20)

pp <- cbind(ind=rind, as.data.frame(alldata[rind, ]$prect)) %>%
  pivot_longer(-ind, values_to="precmd", names_to="day") %>%
  mutate(prect = as.numeric(sub("prec", "", day)),
         precm = 1,
         slope=0,
         faults=0,
         rivers=0,
         lith=alldata$lith[1],
         planc=0,
         profc=0,
         sux=0,
         suy=0,
         timeID=1)

# now making the predictions and extracting the first column
pp$pr <- predict(b_sofr, type="terms", newdata=pp, discrete=FALSE, terms="s(prect):precm")

precm_data <- ggplot(pp) +
  geom_line(aes(x=prect, y=precmd,
                # force discrete colour scheme
                colour=as.factor(ind),
                group=as.factor(ind))) +
  labs(x="Day", y="precm", colour="ind") +
  # give informative labels
  scale_x_continuous(expand=expansion(0, 0)) +
  theme_minimal()
precm_data

precm_data_smooth <- ggplot(pp) +
  geom_line(aes(x=prect, y=pr*precmd,
                colour=as.factor(ind),
                group=as.factor(ind))) +
  labs(x="Day", y="s(day) × precm", colour="ind") +
  # give informative labels
  scale_x_continuous(expand=expansion(0, 0)) +
  theme_minimal()

pp <- pp %>%
  group_by(ind) %>%
  mutate(cpr = cumsum(pr*precmd))

# now build the plot
precm_cumul <- ggplot(pp) +
  geom_line(aes(x=prect, y=cpr,
                colour=as.factor(ind),
                group=as.factor(ind))) +
  labs(x="Day",
       y="Cumulative effect of precm", colour="ind") +
  # give informative labels
  scale_x_continuous(expand=expansion(0, 0)) +
  theme_minimal()


precm_figure <- (precm_data + precm_term)/(precm_data_smooth + precm_cumul) +
  plot_layout(guides="collect")
precm_figure

ggsave(precm_figure, width=8, height=8, file="plots/sofr_precm.pdf")


