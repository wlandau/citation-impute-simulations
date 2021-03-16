
model1 <- function(n, k, density) {

  B <- matrix(0, nrow = k, ncol = k)
  diag(B) <- 0.8
  B[upper.tri(B)] <- 0.3

  directed_dcsbm(
    n = n,
    B = B,
    expected_density = density
  )
}

full_svds <- function(A, k) {
  svds(A, k)
}

zero_imputed_svds <- function(A, k) {
  A_obs <- as(triu(A) * 1, "dgCMatrix")
  svds(A_obs, k)
}

symmetric_svd <- function(A, k) {
  A_obs <- as(triu(A) * 1, "dgCMatrix")
  svds((A_obs + t(A_obs)) / 2, k)
}

cite_impute <- function(A, k) {
  A_obs <- as(triu(A) * 1, "dgCMatrix")
  citation_impute(A_obs, rank = k, max_iter = 20)
}

subspace_loss <- function(svd_true, svd_estimate, params) {
  tibble(
    u_loss = norm(svd_estimate$u - svd_true$u, type = "F") / 2,
    v_loss = norm(svd_estimate$v - svd_true$v, type = "F") / 2
  ) %>%
    bind_cols(params)
}

summarize_losses_model1 <- function(combined_losses) {
  combined_losses %>%
    gather(loss_type, loss, contains("loss")) %>%
    group_by(estimator, n, k, density, loss_type) %>%
    summarize_at(vars(loss), list(~mean(.), ~sd(.))) %>%
    ungroup() %>%
    mutate(
      estimator = as.factor(
        case_when(
          estimator == "loss_full_svds" ~ "Fully Observed SVD",
          estimator == "loss_symmetric_svd" ~ "Symmetrized SVD",
          estimator == "loss_zero_imputed_svds" ~ "Zero Imputed SVD"
        )
      )
    )
}

plot_model1_results <- function(tidied, loss_type = c("u", "v")) {

  loss_type <- match.arg(loss_type)

  label_equals <- function(labs) {
    label_both(labs, sep = " = ")
  }

  plot <- tidied %>%
    filter(loss_type == !!glue("{loss_type}_loss")) %>%
    ggplot() +
    aes(x = n, y = mean, color = estimator, label = estimator) +
    geom_line() +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd)) +
    scale_color_viridis_d(begin = 0.15, end = 0.85, direction = -1) +
    labs(
      title = "Performance of subspace estimators",
      y = glue("{toupper(loss_type)} Loss"),
      x = "Nodes",
      color = "Estimator"
    ) +
    facet_grid(
      rows = vars(k),
      cols = vars(density),
      scales = "free",
      labeller = label_equals
    ) +
    theme_minimal_grid(13) +
    theme(
      panel.spacing = unit(0.7, "lines")
    )

  model1_results_figure_path <- here(
    glue("figures/model1_{loss_type}_loss_results.pdf")
  )

  ggsave(
    model1_results_figure_path,
    plot = plot,
    width = 8.5,
    height = 8.5 * (9/16)
  )

  model1_results_figure_path
}


