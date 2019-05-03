#' Fit multiple models
#'
#' Fits a user-specified model to each imputed dataset.
#'
#' @param imputations A nested tibble, representing multiply imputed datasets.
#' @param model_function A function that takes a dataset as input and returns a fitted
#' model object as output.
#' @return A tibble with the fits of each model.
#' @export
fit_model <- function(imputations, model_function) {
    imps_model <- imputations %>% dplyr::mutate(fit = purrr::map(.data$data, ~my_model(data = .x)),
        tidied = purrr::map(.data$fit, broom::tidy)) %>% tidyr::unnest(.data$tidied)
    imps_model
}

#' Combine Multiple Imputations
#'
#' Pools inferences for regression coefficients for linear models and generalized linear
#' models.
#' @param model_fits An object returned by the function \code{fit_model}
#' @param by_significance Should the p-values be ordered in ascending order?
#' @return A tibble with inferences obtained by multiple imputation combining rules.
#' @export
pool_inferences <- function(model_fits, by_significance = FALSE) {
    pooled_inferences <- model_fits %>% dplyr::rename(est = .data$estimate,
                                                      se = .data$std.error) %>%
        dplyr::group_by(.data$term) %>%
        dplyr::summarise(estimate = mean(.data$est), b = var(.data$est),
                         v = mean(.data$se^2), n = dplyr::n(),
                         std.error = sqrt(.data$v + (.data$n + 1)/.data$n * .data$b),
                         rd = (.data$n + 1)/.data$n * .data$b/.data$v,
                         df = (.data$n - 1) * (1 + 1/.data$rd)^2,
                         statistic = .data$estimate/.data$std.error,
            p.value = 2 * pt(abs(.data$statistic), df = .data$df, lower.tail = F)) %>%
        dplyr::select(.data$term, .data$estimate, .data$std.error, .data$statistic,
                      .data$p.value)
    if (by_significance) {
        pooled_inferences <- pooled_inferences %>% dplyr::arrange(.data$p.value)
    }
    pooled_inferences
}
