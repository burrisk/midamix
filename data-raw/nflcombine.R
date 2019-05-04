library(rvest)
library(dplyr)

year_begin <- 1987
year_end <- 2019

extract_combine_data <- function(url){
  raw <- read_html(url) %>%
    html_nodes('td') %>%
    html_text()
  mat <- matrix(raw, ncol = 13, byrow = T)
  n <- nrow(mat)
  colnames(mat) <- mat[1, ]
  mat <- mat[-c(1, n), ]
}

combine_dfs <- lapply(year_begin:year_end, function(year){
  url <- paste(c("http://nflcombineresults.com/nflcombinedata.php?year=",
                 year, "&pos=&college="), collapse = "")
  extract_combine_data(url)
})

combine_mat <- do.call(rbind, combine_dfs)

to_numeric <- c("Height (in)", "Weight (lbs)", "40 Yard", "Vert Leap (in)",
                "Broad Jump (in)", "Shuttle", "3Cone")
to_integer <- c("Wonderlic", "Bench Press")
to_factor <- c("Name", "College", "POS")
to_ordered <- c("Year")

nfl_combine <- combine_mat %>%
  as_tibble() %>%
  mutate_at(to_numeric, as.numeric) %>%
  mutate_at(to_integer, as.integer) %>%
  mutate_at(to_factor, as.factor) %>%
  mutate_at(to_ordered, as.ordered) %>%
  rename(year = Year, name = Name, college = College, position = POS,
         height = 'Height (in)', weight = "Weight (lbs)", wonderlic = Wonderlic,
         forty_yard_dash = "40 Yard", bench_press = "Bench Press",
         vertical_jump = "Vert Leap (in)",  broad_jump = "Broad Jump (in)",
         shuttle = "Shuttle", three_cone = "3Cone")

nines_to_na <- function(x){
  ifelse(x == 9.99, NA, x)
}

nine_vars <- c("forty_yard_dash", "shuttle", "three_cone")

nfl_combine <- nfl_combine %>%
  mutate_at(nine_vars, nines_to_na)

usethis::use_data(nfl_combine, overwrite = TRUE)


