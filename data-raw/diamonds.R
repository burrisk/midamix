diamonds <- ggplot2::diamonds

set.seed(314)
diamonds <- purrr::map_df(diamonds, function(x) {x[sample(c(TRUE, NA),
                                                          prob = c(0.8, 0.2),
                                                          size = length(x),
                                                          replace = TRUE)]}) %>%
  dplyr::sample_n(1000)

usethis::use_data(diamonds)
