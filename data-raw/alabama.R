library(tidyverse)

alabama <- read_csv("data-raw/bama_raw.csv")

na_to_zero <- function(x){
  ifelse(is.na(x), 0, x)
}

na_to_one <- function(x){
  ifelse(is.na(x), 1, x)
}

zeroed_variables <- c("Income", "HrsWorked")
oned_variables <- c("Educ")

alabama <- alabama %>%
  mutate(Age = as.integer(Age),
         Educ = as.integer(Educ),
         Income = as.integer(Income),
         HrsWorked = as.integer(HrsWorked),
         Medicare = case_when(
           Medicare == 1 ~ 1,
           Medicare == 2 ~ 0
         ),
         Disability = case_when(
           Disability == 1 ~ 1,
           Disability == 2 ~ 0
         ),
         Medicaid = case_when(
           Medicaid == 2 ~ 0,
           Medicaid == 1 ~ 1
         )) %>%
  mutate_at(zeroed_variables, na_to_zero) %>%
  mutate_at(oned_variables, na_to_one) %>%
  dplyr::filter(Income >= 0, Age >= 65 | Medicare == 0 | Disability == 1)

set.seed(314)
validator <- function(y){
  impossible <- (y[5] >= 3 & y[1] <= 3) |
    (y[5] >= 16 & y[1] <= 11) |
    (y[5] >= 21 & y[1] <= 15) |
    (y[5] >= 24 & y[1] <= 19) |
    (y[1] < 65 & y[2] == 1 & y[4] == 0) |
    (y[6] > 100000 & y[3] == 1) |
    (y[7] >= 1 & y[6] <= 0)
  !(impossible)
}

alabama <- alabama[apply(alabama, 1, validator), ] %>%
  dplyr::sample_frac(0.1)

usethis::use_data(alabama)
