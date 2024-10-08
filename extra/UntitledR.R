library(visTimeEvent)
library(survival)
library(dplyr)
a <- survival::pbc
b <- survival::lung

b <- b %>%
  mutate(status = status-1,
         sex1 = ifelse(sex == 1, "male", "female"),
         sex1 = as.factor(sex1),
         status = as.integer(status))
table(b$status, useNA = "always")

result <- km_grouped(data = b,
                     time = "time",
                     event = "status",
                     group = "sex",
                     title = "Test",
                     time_survival = 50,
                     unit = "days",
                     show_label = "median")



library(survival)
library(dplyr)
test_dat <- survival::pbc
test_dat <- test_dat %>%
  mutate(status1 = ifelse(status == 2, 1, 0),
         trt = ifelse(trt == 1, "treatment", "placebo")) %>%
  filter(!is.na(trt))

km_single(data = test_dat,
                     time = "time",
                     event = "status1",
                     title = "Test",
                     time_survival = 50,
                     unit = "days",
                     show_label = "none",
          colors = "blue")
