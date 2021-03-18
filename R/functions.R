flip_coin <- function(p = 0.5) {
  rbinom(1, 1, p)
}

model1 <- function(n, k, expected_degree) {

  B <- matrix(0, nrow = k, ncol = k)

  diag(B) <- 0.8

  directed_dcsbm(
    theta_in = rep(1, n),
    theta_out = rep(1, n),
    B = B,
    expected_density = expected_degree / n
  )
}

model2 <- function(n, k, expected_degree) {

  B <- matrix(0, nrow = k, ncol = k)

  for (i in 1:k) {
    for (j in 1:i) {

      if (flip_coin()) {
        B[i, j] <- 0.4
      } else {
        B[j, i] <- 0.4
      }
    }
  }

  diag(B) <- 0.8

  directed_dcsbm(
    theta_in = rep(1, n),
    theta_out = rep(1, n),
    B = B,
    expected_density = expected_degree / n
  )
}

model3 <- function(n, k, expected_degree) {

  B <- matrix(0, nrow = k, ncol = k)

  for (i in 1:k) {
    for (j in 1:i) {

      if (flip_coin()) {
        B[i, j] <- 0.4
      } else {
        B[j, i] <- 0.4
      }
    }
  }

  diag(B) <- 0.8

  directed_dcsbm(
    theta_in = 1 + rexp(n, 1/5),
    theta_out = 1 + rexp(n, 1/5),
    B = B,
    expected_density = expected_degree / n
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

plot_matrix <- function(X, ...) {
  as_tibble(as.matrix(X), rownames = "row") %>%
    tidyr::gather(col, value, -row) %>%
    mutate(
      row = forcats::fct_inseq(row),
      col = forcats::fct_inseq(col)
    ) %>%
    ggplot(aes(x = col, y = row, fill = value)) +
    geom_raster() +
    scale_fill_gradient2() +
    theme_minimal(16) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank()
    )
}

plot_expectation <- function(model, parameters) {
  E <- model$X %*% tcrossprod(model$S, model$Y)
  plot <- plot_matrix(E)

  n <- parameters$n
  k <- parameters$k
  expected_degree <- parameters$expected_degree

  if (!dir.exists(here("figures/model1_expectations")))
    dir.create(here("figures/model1_expectations"))

  path <- here(
    glue("figures/model1_expectations/n={n}_k={k}_expected_degree={expected_degree}.png")
  )

  ggsave(
    path,
    plot = plot,
    width = 8.5,
    height = 8.5 * (9/16)
  )

  path
}

#' Compute sin-theta distance between subspaces
#'
#' @param u An orthogonal basis for a k-dimensional subspace of
#'   n-dimensional space.
#' @param v An orthogonal basis for a k-dimensional subspace of
#'   n-dimensional space.
#'
#' @return
#' @keywords internal
subspace_loss <- function(u, v) {
  # see [1] Vu and Lei 2013 section 2.3
  # and [2] Rohe, Chatterjee, Yu 2011 Annals of Statistics page 1908

  k <- ncol(u)

  s <- svd(crossprod(u, v))
  k - sum(s$d^2)
}

loss_helper <- function(svd_true, svd_estimate, params) {
  tibble(
    u_loss = subspace_loss(svd_estimate$u, svd_true$u),
    v_loss = subspace_loss(svd_estimate$v, svd_true$v)
  ) %>%
    bind_cols(params)
}

summarize_losses_model1 <- function(combined_losses) {
  combined_losses %>%
    gather(loss_type, loss, contains("loss")) %>%
    group_by(estimator, n, k, expected_degree, loss_type) %>%
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
    aes(x = expected_degree, y = mean, color = estimator, label = estimator) +
    geom_line() +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd)) +
    scale_color_viridis_d(begin = 0.15, end = 0.85, direction = -1) +
    labs(
      title = "Performance of subspace estimators",
      y = glue("{toupper(loss_type)} Loss"),
      x = "Expected degree",
      color = "Estimator"
    ) +
    facet_grid(
      rows = vars(k),
      cols = vars(n),
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


