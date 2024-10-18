library(survival)
library(dplyr)
test_dat <- survival::pbc
test_dat <- test_dat %>%
  mutate(status_death = ifelse(status == 2, 1, 0),
         time = time/365,
         trt = as.factor(ifelse(trt == 1, "treatment", "placebo")),
         trt_num = ifelse(trt == "treatment", 1, 0),
         status = as.factor(status)) %>%
  filter(!is.na(trt))


######### Kaplan-Meier plots
# one group

km_single(data = test_dat,
          time = "time",
          event = "status_death",
          title = "Kaplan-Meier curves for one group",
          unit = "years",
          endpoint = "OS",
          colors = "darkblue",
          show_label = "median")

# two groups

km_grouped(data = test_dat,
          time = "time",
          event = "status_death",
          group = "trt",
          time_survival = 6,
          test = "wilcoxon",
          title = "Kaplan-Meier curves for two groups",
          unit = "years",
          endpoint = "OS",
          colors = c("darkblue", "darkgreen"),
          show_label = "probability")


# weighted, two groups
# CAVE: geht noch nicht: vgl auch MULTIPLY
library(WeightIt)

test_dat_iptw <- test_dat %>%
  filter_at(vars(sex, platelet, trig, age), all_vars(!is.na(.)))

weight_obj <- weightit(trt_num ~ sex + platelet + trig + age, data = test_dat_iptw,
                       method = "glm")

test_dat_iptw$weights <- weight_obj$weights

km_grouped_weighted(data = test_dat_iptw,
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
          x_breaks = c(0, 2, 4, 15))



###############
# competing risk curves
################
# one group
comp_risk_single(data = test_dat,
                 time = "time",
                 event = "status",
                 title = "Competing risk of a single group",
                 unit = "years",
                 colors = c("darkred", "darkblue"),
                 show_label = "probability",
                 time_survival = 8)

# two groups

comp_risk_grouped(data = test_dat,
                 time = "time",
                 event = "status",
                 group = "trt",
                 title = "Competing risk of a single group",
                 unit = "years",
                 colors = c("darkgreen", "deepskyblue3"),
                 show_label = "probability",
                 time_survival = 8)
