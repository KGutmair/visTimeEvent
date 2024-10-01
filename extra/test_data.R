library(survival)
library(dplyr)
test_dat <- survival::pbc
test_dat <- test_dat %>%
  mutate(status1 = ifelse(status == 2, 1, 0),
         trt = ifelse(trt == 1, "treatment", "placebo")) %>%
  filter(!is.na(trt))

data1 <- test_dat
time1 <- "time"
event1 <- "status1"
time_survival <- 900
group <- "trt"
test <- "log-rank"

table(test_dat$status1)

