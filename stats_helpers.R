plot_prior <- function(prior_samples, prior_name = NULL, trans = NULL, hdci_width = c(0.5, 0.89, 0.95), lowerbound = -Inf, cutoff_hdci = 0.999, max_kde_n = 1e6, min_kde_n = 1e3) {

  # Apply transformation to prior samples 
  if (!is.null(trans)) {
    xlab <- paste0(trans, ifelse(is.null(prior_name), "(prior)", paste0("(", prior_name, ")")))
    prior_samples <- do.call(trans, list(prior_samples))
  } else {
    xlab <- ifelse(is.null(prior_name), "(prior)", prior_name)
  }
  
  # Find KDE cutoff value and lenght of x grid (n)
  kde_cutoff <- ggdist::hdci(prior_samples, .width = cutoff_hdci)[2]
  kde_n <- max(min(kde_cutoff, max_kde_n), min_kde_n)
  kde_from <- ifelse(lowerbound == -Inf, floor(min(prior_samples)), lowerbound)
  # Get kernel density estimate (KDE) of prior samples
  dens <- density(prior_samples, n = kde_n, from = kde_from, to = kde_cutoff)

  # Exclude values lower than lowerbound (e.g. KDE will produce positive density at negative x from prior samples close to zero)
  data <- data.frame(x = dens$x, y = dens$y)

  # Get highest density continuous interval for each given width
  # Note: this is not neccessarily the HDI if the prior samples distribution is multi-modal
  data_hdci <- data.frame()
  hdci_width <- sort(hdci_width, decreasing = TRUE) 
  for (i in 1:length(hdci_width)) {
    width_i <- hdci_width[i]
    hdci_i <- ggdist::hdci(prior_samples, .width = width_i)
    data_hdci_i <- data |> 
      dplyr::filter(x >= hdci_i[1] & x <= hdci_i[2]) |> 
      dplyr::mutate(HDCI = width_i, HDCI = as.factor(HDCI))
    data_hdci <- rbind(data_hdci, data_hdci_i)
  }
  # Plot the prior samples distribution and shade HDCI intervals
  ggplot2::ggplot() + 
    ggplot2::geom_area(ggplot2::aes(x, y, fill = HDCI), data = data_hdci, position = "identity") + 
    ggplot2::labs(x = xlab, y = "Density") +
    ggplot2::geom_line(ggplot2::aes(x, y), data = data) + 
    ggplot2::scale_fill_manual(values = ggsci::pal_d3()(length(hdci_width)))
}
plot_priors <- function(brmsfit_priorsonly, plot_params, trans = NULL, lowerbound = -Inf, patch_plots = FALSE, patch_ncol = 2) {

  plots <- list()
  prior_samples <- brmsfit_priorsonly |> as_tibble()

  for (i in 1:length(plot_params)) {
    if (is.null(trans)) {trans_i <- NULL} else {trans_i <- trans[i]}
    if (lowerbound == -Inf) {lowerbound_i <- -Inf} else {lowerbound_i <- lowerbound[i]}
    plots[[i]] <- plot_prior(
      prior_samples = prior_samples |> pull(plot_params[i]), 
      prior_name = plot_params[i], 
      trans = trans_i, 
      lowerbound = lowerbound_i
    )
  }
  # Return either plot list or plots patched together in a single plot object
  if (patch_plots == TRUE) {patchwork::wrap_plots(plots, ncol = patch_ncol)} else {plots}
}