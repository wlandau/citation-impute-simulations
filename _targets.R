library(targets)
library(tarchetypes)
library(magrittr)
library(glue)

options(
  tidyverse.quiet = TRUE,
  clustermq.scheduler = "multicore"
)

# requires dev fastRG
# remotes::install_github("RoheLab/fastRG")

tar_option_set(
  packages = c(
    "fastRG",
    "Matrix",
    "RSpectra",
    "cowplot",
    "glue",
    "tidyverse",
    "here",
    "fastadi"
  )
)

source(here::here("R/functions.R"))

library(logger)
log_threshold(WARN) # hack to hide citation_impute() output

models_chr <- paste0("model", 1:2)

estimators_chr <- c(
  "full_svds",
  # "zero_imputed_svds",
  "symmetric_svd" #,
  # "cite_impute"
)

model_syms <- rlang::syms(models_chr)
estimator_syms <- rlang::syms(estimators_chr)

target_models <- tar_map(

  unlist = FALSE,

  values = tibble::tibble(
    model = model_syms
  ),

  tar_target(
    distribution,
    model(parameters$n, parameters$k, parameters$expected_degree),
    pattern = map(parameters),
    iteration = "list"
  ),

  tar_target(
    realizations,
    purrr::map(1:num_reps, ~sample_sparse(distribution)),
    pattern = map(distribution),
    iteration = "list"
  ),

  tar_target(
    svd,
    svds(distribution),
    pattern = map(population)
  )
)


estimate_values <- tidyr::expand_grid(
  model = models_chr,
  estimator = estimators_chr,
) %>%
  dplyr::mutate(
    estimate = glue("estimates_{model}_{estimator}"),
    data = glue("realizations_{model}")
  ) %>%
  dplyr::mutate_all(rlang::syms)

target_estimates <- tar_eval(
  values = estimate_values,

  tar_target(
    estimate,
    purrr::map(data, ~estimator(.x, k = parameters$k)),
    pattern = map(data)
  )
)

loss_values <- tidyr::expand_grid(
  model = models_chr,
  estimator = estimators_chr,
) %>%
  dplyr::mutate(
    svd = glue("svd_{model}_{estimator}"),
    estimate = glue("estimates_{model}_{estimator}"),
    loss = glue("loss_{estimator}_{model}")
  ) %>%
  dplyr::mutate_all(rlang::syms)

target_loss <- tar_eval(

  values = loss_values,

  tar_target(
    loss,
    purrr::map_dfr(estimate, ~loss_helper(svd, .x, params = parameters)),
    pattern = map(svd, estimate, parameters)
  )
)

target_combined <- tar_combine(
  combined_losses,
  target_estimates[2],
  command = dplyr::bind_rows(!!!.x, .id = "estimator")
)

list(

  tar_target(num_reps, 2),

  tar_target(n, c(500, 1000)), #, 2500)),

  tar_target(expected_degree, c(10, 25, 50, 100, 200)),

  tar_target(k, c(2, 5, 10)),

  tar_target(
    parameters,
    tibble::tibble(n = n, k = k, expected_degree = expected_degree),
    pattern = cross(n, k, expected_degree)
  ),

  # tar_target(
  #   expectation_plots,
  #   plot_expectation(population_graph, parameters),
  #   pattern = map(population_graph, parameters),
  #   format = "file"
  # ),

  target_models,

  target_estimates,

  target_loss
#
#   target_combined,
#
#   tar_target(
#     model1_summarized_losses,
#     summarize_losses_model1(combined_losses)
#   ),
#
#   tar_target(
#     model1_u_loss_plot,
#     plot_model1_results(model1_summarized_losses, loss_type = "u"),
#     format = "file"
#   ),
#
#   tar_target(
#     model1_v_loss_plot,
#     plot_model1_results(model1_summarized_losses, loss_type = "v"),
#     format = "file"
#   )
)

