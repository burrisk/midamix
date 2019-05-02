#' Fits a user-specified model to each imputed dataset
#'
#' @param imputations A nested tibble, representing multiply imputed datasets
#' @param model_function A function that takes a dataset as input and returns a fitted
#' model object as output.
#' @return A tibble with the fits of each model.
fit_model <- function(imputations, model_function) {
    imps_model <- imputations %>% dplyr::mutate(fit = purrr::map(data, ~my_model(data = .x)),
        tidied = purrr::map(fit, broom::tidy)) %>% tidyr::unnest(tidied)
    imps_model
}

#' Pools inferences for regression coefficients for linear models and generalized linear
#' models
#' @param model_fits An object returned by the function \code{fit_model}
#' @param by_significance Should the p-values be ordered in ascending order?
#' @return A tibble with inferences obtained by multiple imputation combining rules.
pool_inferences <- function(model_fits, by_significance = FALSE) {
    pooled_inferences <- model_fits %>% dplyr::rename(est = estimate, se = std.error) %>% dplyr::group_by(term) %>%
        dplyr::summarise(estimate = mean(est), b = var(est), v = mean(se^2), n = n(), std.error = sqrt(v +
            (n + 1)/n * b), rd = (n + 1)/n * b/v, df = (n - 1) * (1 + 1/rd)^2, statistic = estimate/std.error,
            p.value = 2 * pt(abs(statistic), df = df, lower.tail = F)) %>% dplyr::select(term,
        estimate, std.error, statistic, p.value)
    if (by_significance) {
        pooled_inferences <- pooled_inferences %>% dplyr::arrange(p.value)
    }
    pooled_inferences
}
