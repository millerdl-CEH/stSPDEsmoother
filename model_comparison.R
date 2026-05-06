# model comparison

library(ggplot2)
library(mgcv)
library(INLA)
source("smooth.construct.spdeST.smooth.spec.R")
library(gratia)
library(dplyr)
library(tidyr)

# load the fitted distributed lag model
load("dl_model.rda")
# SoFR model
load("sofr_model.rda")
load("data.rda")

# SoFR model seems to be better by AIC at least?
AIC(b_dl, b_sofr)

# not much difference between models, so interpretation is important?
preds <- data.frame(dl=predict(b_dl, type="response"),
                    sofr=predict(b_sofr, type="response"))

ggplot(preds, aes(x=dl, y=sofr)) +
  geom_point() +
#  geom_density_2d() +
  labs(x="Distributed lag", y="Scalar-on-function regression") +
  theme_minimal()


# compare parameter estimates for non-DL/SoFR parameters
coef_dl <- coef(b_dl)
coef_sofr <- coef(b_sofr)

coef_dl <- coef_dl[!grepl("precm", names(coef_dl))]
coef_sofr <- coef_sofr[!grepl("precm", names(coef_sofr))]

plot(coef_dl, coef_sofr, asp=1)
abline(a=0, b=1, col="red")

## compare effects

sm <- rbind(
  smooth_estimates(b_dl, select=c("s(slope)", "s(faults)", "s(rivers)")) |>
    add_confint() |>
    mutate(model="dl"),
  smooth_estimates(b_sofr, select=c("s(slope)", "s(faults)", "s(rivers)")) |>
    add_confint() |>
    mutate(model="sofr")) |>
  select(.smooth, .estimate, slope, rivers, faults, .lower_ci, .upper_ci, model) |>
  pivot_longer(cols=c(slope, rivers, faults), names_to="covariate", values_to="x")

p_eff <- ggplot(sm) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci, x = x,
                  group=model, fill=model), alpha = 0.2) +
  geom_line(aes(x = x, y = .estimate, group=model, colour=model), lwd = 1.2) +
  facet_wrap(~.smooth) +
  labs(y = "Partial effect", fill="Model", colour="Model", x="Covariate value") +
  theme_minimal()

p_eff
ggsave(p_eff, file="plots/compare_smooth_effects.pdf", width=10, height=12)


aa <- smooth_estimates(b_dl, select=c("prect"), partial_match=TRUE) |>
  add_confint()
bb <- smooth_estimates(b_sofr, select=c("prect"), partial_match=TRUE) |>
  add_confint()


