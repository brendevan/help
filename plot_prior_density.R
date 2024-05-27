# ===================================================
#             PLOT PRIOR SIMULATIONS
# ===================================================
# Plot density of n_samples samples from the given prior distribution and shade the given highest continuous density intervals (HDCIs). To add a lowerbound to the prior samples (e.g. as per the default brms random effect sd prior with a lowerbound at zero) set truncate_lower. If a link funciton other than the identity link is set, backtransformed prior distributions are also plotted (i.e. on the response scale). There are 2 plot styles defined: (1) plot_prior_density plots a continuous density and (2) plot_prior_density_dots a quantile dotplot of the prior sample density.  

plot_prior_density <- function(
  prior_def, # Prior distribution definition used to generate samples
  truncate_lower = NULL,
  link = "identity",    # Link function for paramater; if supplied, samples are backtransformed using the inverse link function
  type = c("density", "dots"), # Plot type
  interval_widths = c(0.5, 0.85, 0.99),   # Interval widths to plot (HDCI for 'density'; defined by point_interval for 'dots')
  interval_colors = c("#9ECAE1", "#3182BD", "#065085"), 
  cutoff_hdci = 0.999,  # type == density: Density is not estimated for prior values outside the cutoff_hdci HDCI e.g. 0.999 HDCI
  max_kde_n = 1e6,      # type == density: Max number of grid cells at which to estimate the density
  min_kde_n = 1e3,      # type == density: Min number of grid cells at which to estimate the density
  hdci_table = "none",  # type == density: Whether or not to print the table of HDCIs: "none", "print", "return"
  ndots = 100,          # type == dots: Number of dots to plot (sets ggdist::stat_dotsinterval quantiles argument)
  point_interval = "median_hdci", # type == dots: 
  na_color = "#676769", # type == dots: 
  nsamples = 1e5,       # Number of samples over which density and HDCIs are calculated
  patched = TRUE,       # Whether to patch plots together when multiple are produced
  subtitle_width = 80,  # Number of characters after which the subtitle is wrapped to the next line
  adjust_margin = FALSE, # Try to fix the top margin of the plot if subtitle is not fully shown, 
  seed = "random",       # Integer; if "random" (default) a random integer in [0, 1000] is used
  ...                    # Arguments for ggdist functions e.g. alpha = 0.5
) {
  require(ggplot2)
  require(ggdist)
  require(ggtext)
  require(latex2exp)
  require(dplyr)
  require(stringr)
  require(patchwork)
  require(rlang)
  type <- rlang::arg_match(type)

  # Parse parameter name from prior def
  param <- .get_param_name(prior_def)
  
  # Get samples from prior distribution  
  if (seed == "random") seed <- round(runif(1)*1000, 0)
  message(paste("Random seed =", seed))
  prior_samples <- .get_prior_samples(prior_def, nsamples, seed, lower = truncate_lower, msg = TRUE)

  # Plot samples
  plots <- list()
  title <- ifelse(
    param == "sigma", 
    prior_def |> str_split_1("~") |> str_replace("sd", "sigma") |> paste(collapse = "~"), 
    prior_def
  )
  if (type == "density") {
    p1 <- .plot_samples_density(prior_samples, interval_widths, cutoff_hdci, max_kde_n, min_kde_n, hdci_table, subtitle_width = subtitle_width)
  } 
  if (type == "dots") {
    p1 <- .plot_samples_dots(prior_samples, ndots, point_interval, interval_widths, interval_colors, na_color, subtitle_width = subtitle_width)
  }
  plots$prior <- p1 + labs(x = latex2exp::TeX(paste0("\\", param)), title = latex2exp::TeX(paste0("\\", title))
  )

  # Plot Delta ~ N(0, sd) if parameter is random effects sd
  if (param == "sigma") {
    if (any(prior_samples < 0)) stop(
      "negative prior samples produced for random effect sd. 
       Try choosing a 0+ distribution such as a half-distribution or set truncate_lower = 0 to only plot samples >= 0"
    )
    prior_samples <- rnorm(nsamples, mean = 0, sd = prior_samples)
    if (type == "density") {
      p2 <- .plot_samples_density(prior_samples, interval_widths, cutoff_hdci, max_kde_n, min_kde_n, hdci_table, subtitle_width = NA)
    } 
    if (type == "dots") {
      p2 <- .plot_samples_dots(prior_samples, ndots, point_interval, interval_widths, interval_colors, na_color, subtitle_width = NA) 
    }
    plots$random_intercept <- p2 + labs(x = latex2exp::TeX("\\Delta"), title =  latex2exp::TeX("\\Delta ~ Normal(0, \\sigma)"))
    message("sd prior converted to random intercept as RI ~ N(0, sd)")
  }

  # Plot backtransformed samples if not identity-link
  if (link != "identity") {
    inv_link <- .inverse_link(link)
    prior_samples <- .backtransform(prior_samples, inv_link)
    xlab <- if (param == "sigma") paste0(inv_link, "(\\Delta)") else xlab <- paste0(inv_link, "(\\", param, ")")
    if (type == "density") {
      p3 <- .plot_samples_density(prior_samples, interval_widths, cutoff_hdci, max_kde_n, min_kde_n, hdci_table, subtitle_width = NA)
    } 
    if (type == "dots") {
      p3 <- .plot_samples_dots(prior_samples, ndots, point_interval, interval_widths, interval_colors, na_color, subtitle_width = NA) 
    }
    plots$backtransformed <- p3 + 
    labs(x = latex2exp::TeX(xlab)) 
    message(paste0("inverse link function applied [link: ", link, "; inverse link: ", inv_link), "]")
  }
  # =====> RETURN PLOTS
  if (patched) {
    plots <- plots |> 
      patchwork::wrap_plots(ncol = 1)
  }
  browser()
  return(plots)
}


# ==================================
#             HELPERS
# ==================================
.supported_prior_distributions <- function() {
  list(
    `stats::rnorm` = c("normal", "gaussian", "N"),
    `extraDistr::rhnorm` = c("half_normal", "half_gaussian", "half_N"),
    `stats::rlnorm` = c("lognormal", "log_normal", "log_N"),
    `stats::rgamma` = c("gamma"),
    `stats::rexp` = c("exponential", "exp"), 
    `stats::runif` = c("uniform", "U", "unif"), 
    `stats::rpois` = c("poisson"), 
    `brms::rstudent_t` = c("student_t", "t"), 
    `extraDistr::rinvgamma` = c("inv_gamma", "invgamma", "inverse_gamma", "inversegamma")
  )
}


.plot_samples_density <- function(prior_samples, hdci_width, cutoff_hdci, max_kde_n, min_kde_n, hdci_table, subtitle_width) {

  # Exclude HDCI widths greater than the cutoff
  hdci_width <- hdci_width[hdci_width <= cutoff_hdci]
  # Get KDE to and from arguments
  kde_cutoff <- ggdist::hdci(prior_samples, .width = cutoff_hdci)[2]
  kde_from <- floor(ggdist::hdci(prior_samples, .width = cutoff_hdci)[1])
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
  if (hdci_table != "none") {
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
  subtitle_bounds <- if (cutoff_hdci >= 1) "" else paste0("The density is estimated only for values within the ", cutoff_hdci*100, "% HDCI.")
  sub <- paste0(
    "Density plot of ", formatC(length(prior_samples), format="d", big.mark=","), " samples along with Highest Continuous Density Intervals (HDCIs).", 
    subtitle_bounds
  )
  p <- p + 
    ggplot2::labs(
      subtitle = if(!(is.na(subtitle_width) || is.null(subtitle_width))) stringr::str_wrap(sub, subtitle_width) else waiver(),
      x = "Parameter (change this label)", 
      y = "Density"
    ) + 
    ggplot2::theme(legend.text = ggtext::element_markdown())

  return(p)
}

.plot_samples_dots <- function(samples, ndots, point_interval, interval_widths, interval_colors, na_color, subtitle_width) {
  # NOTE: ggdist seems to be unreliable in its point calculations, e.g. 
  # set.seed(710); mode_hdci(exp(rnorm(1e5, 0, rexp(1e5, 2)))
  if (str_detect(point_interval, "mode")) warning("ggdist mode calculations are sometimes suspect and or fail to calculate altogether; beware.")
  # Try calling the point_interval function on the samples to check for errors
  try_point_interval <- tryCatch(
    do.call(point_interval, list(samples)),
    error = function(e) NULL
  )
  if (is.null(try_point_interval)) {
    point_interval = "median_hdci"
    warning(paste0("there was an error calculating", point_interval, "; using median_hdci instead.")) 
  }
  # PLOT STYLING
  gg_remove_yaxis <- theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
  sub <- paste("Quantile dotplot calculated on", formatC(length(samples), format="d", big.mark=","), "samples with dots shaded by Highest Continuous Density Interval (HDCI)")
  # CREATE PLOT
  p <- data.frame(x = samples) |> ggplot()
  p <- p +
    stat_dots(
      aes(x, 
        slab_fill = after_stat(level), 
        slab_color = after_stat(level), 
      ), 
      quantiles = ndots, 
      point_interval = point_interval, 
      .width = interval_widths, 
      layout = "bin"
    ) + 
    scale_fill_manual(
      values = interval_colors, na.value = na_color, aesthetics = c("slab_fill", "slab_color")
    ) +
    # scale_color_manual(
    #   values = interval_colors, na.value = na_color, aesthetics = "slab_color"
    # ) +
    labs(
      x = paste0("Prior samples (n =", ndots, ")"),
      y = "Density",
      slab_fill = "HDCI", slab_color = "HDCI",
      subtitle = if(!(is.na(subtitle_width) || is.null(subtitle_width))) stringr::str_wrap(sub, subtitle_width) else waiver()
    ) + 
    gg_remove_yaxis

  return(p)
}

.get_prior_samples <- function(prior_def, n_samples, seed, lower = NULL, msg = FALSE) {
  dists <- .supported_prior_distributions()
  prior_def <- str_split_1(prior_def, "~")[2] |> str_trim() # e.g. beta ~ N(0, 1) -> N(0, 1)
  distr <- prior_def |> 
    str_extract("^\\w+") |>  # e.g. N(0, 1) -> N
    str_to_lower()
  rfunc <- names(which(sapply(dists, function(x) distr %in% tolower(x)))) # e.g. N -> rnorm
  if (nchar(rfunc) == 0) stop(paste0("prior_def must be one of [", paste(unlist(unlist(dists)), collapse = ", "), "]. Add distributions in .supported_prior_distributions()"))
  args <- prior_def |> str_to_lower() |> 
    str_remove(distr) |> str_remove("\\(") |> str_remove("\\)") |> 
    str_split_1(",") |> str_trim()
  args <- c(list(n = n_samples), as.numeric(as.list(args)))
  rfunc_call <- paste0(rfunc, "(", paste(args, collapse = ", "), ")")
  set.seed(seed)
  prior_samples <- tryCatch(
    eval(parse(text = rfunc_call)), 
    error = function(e) stop(paste0("Trying to call `", rfunc_call, "`. \n", e))
  )
  if (!is.null(lower)) {
    prior_samples <- prior_samples[prior_samples >= lower]
    while (length(prior_samples) < n_samples) {
      prior_samples <- c(prior_samples, tryCatch(
        eval(parse(text = rfunc_call)), 
        error = function(e) stop(paste0("Trying to call `", rfunc_call, "`. \n", e))
      ))
      prior_samples <- prior_samples[prior_samples >= lower]
    }
    prior_samples <- prior_samples[1:n_samples]
  }
  if (!(length(prior_samples) == n_samples)) stop(paste0("Trying to call `", rfunc_call, "` does not generate the expected number of samples."))
  if (msg) message(paste0(length(prior_samples), " prior samples produced via `", rfunc_call, "`", if (!is.null(lower)) paste(" truncated at", lower)))
  return(prior_samples)
}
.get_param_name <- function(prior_def) {
  if (!str_detect(prior_def, "~")) stop("prior_def should be of form '<parameter> ~ <distr_name>(<distr_parameters>)'")
  param <- str_trim(str_split_1(prior_def, "~")[1])
  p <- tolower(param)
  if (p %in% c("sd", "sigma")) param <- "sigma"
  if (param %in% c("b", "beta")) param <- "beta"
  return(param)
}
.inverse_link <- function(link) {
  link <- tolower(link)
  inv_link <- NULL
  if (link == "identity") inv_link <- "identity"
  if (link == "log") inv_link <- "exp"
  if (is.null(inv_link)) stop("link function not supported; implement in .inverse_link()")
  return(inv_link)
}
.backtransform <- function(samples, inverse_link) {
  # Backtransform samples back to the response scale=
  eval(parse(text = paste0(inverse_link, "(samples)")))
}