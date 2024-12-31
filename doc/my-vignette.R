## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(visTimeEvent)
library(survival)
library(dplyr)

## ----example------------------------------------------------------------------
test_dat <- survival::pbc
test_dat <- test_dat %>%
  mutate(status_death = ifelse(status == 2, 1, 0),
         status_death = as.integer(status_death),
         time = time/365,
         trt = as.factor(ifelse(trt == 1, "treatment", "placebo")),
         trt_num = ifelse(trt == "treatment", 1, 0),
         status = as.factor(status)) %>%
  filter(!is.na(trt))

## ----KM_plot_one_group, fig.width=7-------------------------------------------
plot_km <- km_single(
  data = test_dat,
  time = "time",
  event = "status_death",
  title = "Kaplan-Meier curves for one group",
  unit = "years",
  endpoint = "OS",
  colors = "darkblue",
  show_label = "median"
)
plot_km[[1]]

## ----KM_plot_two_groups, fig.width=7, fig.height= 5---------------------------
plot_km1 <- km_grouped(
  data = test_dat,
  time = "time",
  event = "status_death",
  group = "trt",
  time_survival = 6,
  #test = "wilcoxon",
  title = "Kaplan-Meier curves for two groups",
  unit = "years",
  endpoint = "OS",
  colors = c("darkblue", "darkgreen"),
  show_label = "probability"
)

plot_km1[[1]]

## ----KM_plot_IPTW_two_groups, fig.width=7, fig.height= 6----------------------
library(WeightIt)

test_dat_iptw <- test_dat %>%
  filter_at(vars(sex, platelet, trig, age), all_vars(!is.na(.)))

weight_obj <- weightit(trt_num ~ sex + platelet + trig + age,
  data = test_dat_iptw,
  method = "glm"
)

test_dat_iptw$weights <- weight_obj$weights

test_dat_iptw <- test_dat_iptw %>%
  mutate(trt = as.factor(if_else(trt == "placebo", 0, 1)))

plot_km_iptw <- km_grouped_weighted(
  data = test_dat_iptw,
  time = "time",
  event = "status_death",
  group = "trt",
  weight_col = "weights",
  group_name = c("treatment", "placebo"),
  endpoint = "OS",
  title = "Weighted KM- curves with two groups",
  unit = "years",
  show_label = "probability",
  time_survival = 6,
  x_lim = c(0, 14),
  x_breaks = c(0, 2, 4, 6, 8, 12, 14),
  risk_table = TRUE,
  color_curves = c("darkgreen", "deepskyblue3"),
  legend_placement = c(0.6, 0.2)
)

plot_km_iptw[[1]]

## ----cumulative_incidence_plot_one_group, fig.width=7-------------------------
plot_cum <- comp_risk_single(
  data = test_dat,
  time = "time",
  event = "status",
  title = "Competing risk of a single group",
  unit = "years",
  colors = c("darkred", "darkblue"),
  show_label = "probability",
  time_survival = 8
)

plot_cum[[1]]

## ----cumulative_incidence_plot_two_groups, fig.width=7, fig.height=5----------
plot_cum1 <- comp_risk_grouped(
  data = test_dat,
  time = "time",
  event = "status",
  group = "trt",
  title = "Competing risk of two group",
  unit = "years",
  colors = c("darkgreen", "deepskyblue3"),
  show_label = "probability",
  time_survival = 8
)
plot_cum1[[1]]

