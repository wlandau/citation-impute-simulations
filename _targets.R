library(targets)
library(tarchetypes)

options(
  tidyverse.quiet = TRUE,
  clustermq.scheduler = "multicore"
)

# requires dev fastRG
# remotes::install_github("RoheLab/fastRG")

tar_option_set(
  packages = c("fastRG", "Matrix", "RSpectra", "cowplot", "glue", "tidyverse", "here")
)

source(here::here("R/functions.R"))

##### TODO list ----------------------------------------------------------------
#
# - Generalize to more than one generative model
#
##### END TODO list ------------------------------------------------------------

target_estimates <- tar_map(

  unlist = FALSE,

  values = tibble::tibble(
    estimator = rlang::syms(
      c(
        "full_svds",
        "zero_imputed_svds",
        "symmetric_svd"
      )
    )
  ),

  tar_target(
    estimate,
    purrr::map(model1_realizations, ~estimator(.x, k = model1_parameters$k)),
    pattern = map(model1_realizations, model1_parameters)
  ),

  tar_target(
    loss,
    purrr::map_dfr(estimate, ~subspace_loss(population_svd, .x, params = model1_parameters)),
    pattern = map(population_svd, estimate, model1_parameters)
  )
)

target_combined <- tar_combine(
  combined_losses,
  target_estimates[2],
  command = dplyr::bind_rows(!!!.x, .id = "estimator")
)

num_reps <- 10

list(

  tar_target(n, c(100, 250, 500, 1000, 2500)),

  tar_target(k, c(2, 5, 10)),

  tar_target(density, c(0.1, 0.2, 0.3, 0.4)),

  tar_target(
    model1_parameters,
    tibble::tibble(n = n, k = k, density = density),
    pattern = cross(n, k , density)
    # tidyr::expand_grid(n, k, density)
  ),

  tar_target(
    population_graph,
    model1(model1_parameters$n, model1_parameters$k, model1_parameters$density),
    pattern = map(model1_parameters),
    iteration = "list"
  ),

  tar_target(
    model1_realizations,
    purrr::map(1:num_reps, ~sample_sparse(population_graph)),
    pattern = map(population_graph),
    iteration = "list"
  ),

  tar_target(
    population_svd,
    svds(population_graph),
    pattern = map(population_graph)
  ),

  target_estimates,

  target_combined,

  tar_target(
    model1_summarized_losses,
    summarize_losses_model1(combined_losses)
  ),

  tar_target(
    model1_u_loss_plot,
    plot_model1_results(model1_summarized_losses, loss_type = "u")
  ),

  tar_target(
    model1_v_loss_plot,
    plot_model1_results(model1_summarized_losses, loss_type = "v")
  )
)

