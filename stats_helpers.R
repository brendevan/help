plot_prior_density <- function(
  parameter_name,      # The name of the model parameter for which the prior pertains; used for plot title and labs
  prior_distribution,  # the prior distribution e.g. "normal"
  prior_args,          # arguments for the distribution e.g. list(mean = 0, sd = 1)
  n_samples = 1e6,     # how many samples to draw from prior distribution to form density
  trans = NULL,        # transformation to be applied to samples before plotting
  hdci_width = c(0.5, 0.89, 0.95, 0.99),   # HDCI widths to plot
  lowerbound = NULL,   # Lowerbound for density estimates, if not supplied cutoff_hdci is used
  cutoff_hdci = 0.999, # Density is not estimated for prior values outside the cutoff_hdci HDCI e.g. 0.999 HDCI
  max_kde_n = 1e6,     # Max number of grid cells at which to estimate the density
  min_kde_n = 1e3,     # Min number of grid cells at which to estimate the density
  subtitle_width = 60,  # Number of characters after which the subtitle is wrapped to the next line
  adjust_margin = FALSE # Try to fix the top margin of the plot if subtitle is not fully shown
) {

  # Set values based on prior_distribution
  dists <- list(
    normal = "rnorm", 
    gamma = "rgamma"
  )
  arg_strings <- c()
  for (i in 1:length(prior_args)) {
    name <- names(prior_args)[i]
    value <- prior_args[[i]]
    arg_strings[i] <- paste(name, "=", value)
  }
  dist_string <- paste0(parameter_name, " ~ ", stringr::str_to_title(prior_distribution), "(", paste(arg_strings, collapse = ", "), ")")

  # Generate prior samples
  if (prior_distribution %in% names(dists)) {
    dist_func <- dists[[prior_distribution]]
  } else {
    stop(paste0("Unsupported prior_dsitribution. Supported distributions are: ", paste(names(dists), collapse = " ")))
  }
  prior_samples <- do.call(dist_func, args = append(list(n_samples), prior_args))

  # Apply transformation to prior samples 
  if (!is.null(trans)) {
    xlab <- paste0(trans, "(", parameter_name, ")")
    prior_samples <- do.call(trans, list(prior_samples))
  } else {
    xlab <- parameter_name
  }
  
  # Exclude HDCI widths greater than the cutoff
  hdci_width <- hdci_width[hdci_width <= cutoff_hdci]

  # Find KDE cutoff value and lenght of x grid (n)
  kde_cutoff <- ggdist::hdci(prior_samples, .width = cutoff_hdci)[2]
  kde_n <- max(min(kde_cutoff, max_kde_n), min_kde_n)
  kde_from <- ifelse(is.null(lowerbound), 
    floor(ggdist::hdci(prior_samples, .width = cutoff_hdci)[1]), 
    lowerbound
  )
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

  # Message HDCI intervals
  hdci_tbl <- data.frame()
  for (i in 1:length(hdci_width)) {
    width_i <- hdci_width[i]
    hdci_i <- ggdist::hdci(prior_samples, .width = width_i)
    hdci_tbl_i <- data.frame(`HDCI width` = width_i, From = hdci_i[1], To = hdci_i[2])
    hdci_tbl <- rbind(hdci_tbl, hdci_tbl_i)
  }
  # Plot the prior samples distribution and shade HDCI intervals
  p <- ggplot2::ggplot() + 
    ggplot2::geom_area(ggplot2::aes(x, y, fill = HDCI), data = data_hdci, position = "identity") + 
    ggplot2::geom_line(ggplot2::aes(x, y), data = data) + 
    ggplot2::scale_fill_manual(values = ggsci::pal_d3()(length(hdci_width)))

  # Add subtitle to plot
  if (!is.null(trans)) {
    trans_phrase <- paste0(trans, "-transformed ")
    if (trans == "exp") trans_phrase <- "exponentiated "
  } else {
    trans_phrase <- ""
  }
  if (cutoff_hdci >= 1 & is.null(lowerbound)) {
    subtitle_bounds <- ""
  } else if (cutoff_hdci >= 1 & !is.null(lowerbound)) {
    subtitle_bounds <- paste0("The density is estimated only for prior values greater than the lowerbound of ", lowerbound, ".")
  } else if (cutoff_hdci < 1 & is.null(lowerbound)) {
    subtitle_bounds <- paste0("The density is estimated only for prior values within the ", cutoff_hdci*100, "% HDCI.")
  } else if (cutoff_hdci < 1 & !is.null(lowerbound)) {
    subtitle_bounds <- paste0("The density is estimated only for prior values within the ", cutoff_hdci*100, "% HDCI and greater than the lowerbound of ", lowerbound, ".")
  }
  subtitle <- paste0(
    "Density plot of ", formatC(length(prior_samples), format="d", big.mark=","), " ",
    trans_phrase, "samples of ", parameter_name, " from the distribution ", dist_string,  
    " along with Highest Continuous Density Intervals (HDCIs).", 
    subtitle_bounds
  )
  p <- p + ggplot2::labs(
    subtitle = stringr::str_wrap(subtitle, subtitle_width),
    x = xlab, 
    y = "Density"
  )

  # Adjust plot margins to account for wrapped subtitle
  if (adjust_margin) {
    subtitle_nlines <- ceiling(nchar(subtitle)/subtitle_width)
    p <- p + ggplot2::theme(plot.margin = grid::unit(c(subtitle_nlines, 0, 0, 0), "lines"))
  }
  
  # Print HDCIs to console as dataframe
  print(hdci_tbl)

  return(p)
}

plot_prior_density2 <- function(
  prior_samples,
  hdci_width = c(0.5, 0.89, 0.95, 0.99),   # HDCI widths to plot
  lowerbound = NULL,   # Lowerbound for density estimates, if not supplied cutoff_hdci is used
  cutoff_hdci = 0.999, # Density is not estimated for prior values outside the cutoff_hdci HDCI e.g. 0.999 HDCI
  max_kde_n = 1e6,     # Max number of grid cells at which to estimate the density
  min_kde_n = 1e3,     # Min number of grid cells at which to estimate the density
  subtitle_width = 60, # Number of characters after which the subtitle is wrapped to the next line
  adjust_margin = FALSE, # Try to fix the top margin of the plot if subtitle is not fully shown, 
  hdci_table = "no"    # Whether or not to print the table of HDCIs: "none", "print", "return"
) {
  
  # Exclude HDCI widths greater than the cutoff
  hdci_width <- hdci_width[hdci_width <= cutoff_hdci]

  # Get KDE to and from arguments as:
  #    - from: if lowerbound is given, use lower bound; otherwise use cutoff_hdci
  #    - to:   use cutoff_hdci
  kde_cutoff <- ggdist::hdci(prior_samples, .width = cutoff_hdci)[2]
  kde_from <- ifelse(is.null(lowerbound), 
    floor(ggdist::hdci(prior_samples, .width = cutoff_hdci)[1]), 
    lowerbound
  )
  # Get number of grid cells at which to estimate density
  kde_n <- max(min(ceiling(kde_cutoff), max_kde_n), min_kde_n)
  # Get kernel density estimate (KDE) of prior samples
  dens <- density(prior_samples, n = kde_n, from = kde_from, to = kde_cutoff)
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

  # # HDCIs table (return or print if specified)
  if (hdci_table != "no") {
    hdci_tbl <- data.frame()
    for (i in 1:length(hdci_width)) {
      width_i <- hdci_width[i]
      hdci_i <- ggdist::hdci(prior_samples, .width = width_i)
      hdci_tbl_i <- data.frame(`HDCI width` = width_i, From = hdci_i[1], To = hdci_i[2])
      hdci_tbl <- rbind(hdci_tbl, hdci_tbl_i)
    }
    if (hdci_table == "print") print(hdci_tbl)
    if (hdci_table == "return") return(hdci_table)
  }
  
  # Plot the prior samples distribution and shade HDCI intervals
  p <- ggplot2::ggplot() + 
    ggplot2::geom_area(ggplot2::aes(x, y, fill = HDCI), data = data_hdci, position = "identity") + 
    ggplot2::geom_line(ggplot2::aes(x, y), data = data) + 
    ggplot2::scale_fill_manual(values = ggsci::pal_d3()(length(hdci_width)))

  # Add subtitle to plot
  if (cutoff_hdci >= 1 & is.null(lowerbound)) {
    subtitle_bounds <- ""
  } else if (cutoff_hdci >= 1 & !is.null(lowerbound)) {
    subtitle_bounds <- paste0("The density is estimated only for values greater than the lowerbound of ", lowerbound, ".")
  } else if (cutoff_hdci < 1 & is.null(lowerbound)) {
    subtitle_bounds <- paste0("The density is estimated only for values within the ", cutoff_hdci*100, "% HDCI.")
  } else if (cutoff_hdci < 1 & !is.null(lowerbound)) {
    subtitle_bounds <- paste0("The density is estimated only for values within the ", cutoff_hdci*100, "% HDCI and greater than the lowerbound of ", lowerbound, ".")
  }
  subtitle <- paste0(
    "Density plot of ", formatC(length(prior_samples), format="d", big.mark=","), " samples along with Highest Continuous Density Intervals (HDCIs).", 
    subtitle_bounds
  )
  p <- p + ggplot2::labs(
    subtitle = stringr::str_wrap(subtitle, subtitle_width),
    x = "Parameter (change this label)", 
    y = "Density"
  )

  return(p)
}