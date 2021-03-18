library(targets)
library(tarchetypes)

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

target_models <- tar_map(

  unlist = FALSE,

  values = tibble::tibble(
    model = rlang::syms(paste0("model", 1:2))
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


target_estimates <- tar_map(

  unlist = FALSE,

  values = tibble::tibble(

    estimator = rlang::syms(
      c(
        "full_svds",
        # "zero_imputed_svds",
        "symmetric_svd" #,
        # "cite_impute"
      )
    )
  ),


tar_target(
  estimate,
  purrr::map(realizations, ~estimator(.x, k = parameters$k)),
  pattern = map(realizations, parameters)
)

  # tar_target(
  #   loss,
  #   purrr::map_dfr(estimate, ~loss_helper(population_svd, .x, params = parameters)),
  #   pattern = map(population_svd, estimate, parameters)
  # )
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

  target_estimates
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

