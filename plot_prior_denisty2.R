
plot_prior_density2 <- function(
  prior_samples,
  hdci_width = c(0.5, 0.89, 0.95, 0.99),   # HDCI widths to plot
  lowerbound = NULL,   # Lowerbound for density estimates, if not supplied cutoff_hdci is used
  cutoff_hdci = 0.999, # Density is not estimated for prior values outside the cutoff_hdci HDCI e.g. 0.999 HDCI
  max_kde_n = 1e6,     # Max number of grid cells at which to estimate the density
  min_kde_n = 1e3,     # Min number of grid cells at which to estimate the density
  subtitle_width = 60, # Number of characters after which the subtitle is wrapped to the next line
  adjust_margin = FALSE, # Try to fix the top margin of the plot if subtitle is not fully shown, 
  hdci_table = "none"    # Whether or not to print the table of HDCIs: "none", "print", "return"
) {

  require(ggplot2)
  require(dplyr)
  require(ggdist)
  require(ggtext)
  require(stringr)

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
  digits <- function(x) nchar(as.character(x))
  myround <- function(x) {
    ifelse(x < 1 & x > -1, 
      ifelse(digits(signif(x, 1)) - 2 > 5,  # if number of DPs greater than 5
        round(x, 5),  # use 0.00000
        signif(x, 1)  # otherwise express to 1 sigfig
      ), 
      round(x, 1)     # if not in (-1, 1) round to 1 DP
    )
  }
  data_hdci <- data.frame()
  hdci_width <- sort(hdci_width, decreasing = TRUE) 
  for (i in 1:length(hdci_width)) {
    width_i <- hdci_width[i]
    hdci_i <- ggdist::hdci(prior_samples, .width = width_i)
    data_hdci_i <- data |> 
      dplyr::filter(x >= hdci_i[1] & x <= hdci_i[2]) |> 
      dplyr::mutate(
        HDCI = paste0("**", width_i,"** [", myround(hdci_i[1]), ", ", myround(hdci_i[2]), "]"), 
        HDCI = as.factor(HDCI)
      )
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
  p <- p + 
    ggplot2::labs(
      subtitle = stringr::str_wrap(subtitle, subtitle_width),
      x = "Parameter (change this label)", 
      y = "Density"
    ) + 
    ggplot2::theme(legend.text = ggtext::element_markdown())

  return(p)
}