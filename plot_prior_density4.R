# ==================================
#            DEV NOTES
# ==================================
# - Remove y-axis (label and breaks, text) all together?


# ==================================
#          MAIN FUNCTION
# ==================================
plot_prior_density4 <- function(
  prior_def, # Prior distribution definition used to generate samples
  link = "identity",  # Link function for paramater; if supplied, samples are backtransformed using the inverse link function
  ndots = 100, # Number of dots to plot (sets ggdist::stat_dotsinterval quantiles argument)
  point_interval = "median_hdci",
  interval_widths = c(0.66, 0.97),
  interval_colors = c("#9ECAE1", "#3182BD"), 
  na_color = "#676769",
  vline_color = "#c5b3cf",
  nsamples = 1e5, # Number of samples over which point_interval is calculated
  patched = TRUE, # Whether to patch plots together when multiple are produced
  subtitle_width = 60, # Number of characters after which the subtitle is wrapped to the next line
  adjust_margin = FALSE, # Try to fix the top margin of the plot if subtitle is not fully shown, 
  seed = "random",       # Integer; if "random" (default) a random integer in [0, 1000] is used
  ...                    # Arguments for ggdist::geom_dotsinterval(); e.g. alpha = 0.5
) {
  require(ggplot2)
  require(ggdist)
  require(ggtext)
  require(latex2exp)
  require(dplyr)
  require(stringr)
  require(patchwork)

  param <- .get_param_name(prior_def)
  plots <- list()

  # =====> GET PRIOR SAMPLES AND HDCI
  if (seed == "random") seed <- round(runif(1)*1000, 0)
  message(paste("Random seed =", seed))
  prior_samples <- .get_prior_samples(prior_def, nsamples, seed, msg = TRUE)

  # =====> PLOT PRIOR SAMPLES
  title <- ifelse(param == "sigma", 
    prior_def |> str_split_1("~") |> str_replace("sd", "sigma") |> paste(collapse = "~"), 
    prior_def
  )
  plots$prior <- .plot_samples(prior_samples, ndots, point_interval, interval_widths, interval_colors, na_color) + 
    labs(
      x = latex2exp::TeX(paste0("\\", param)), 
      title = latex2exp::TeX(paste0("\\", title))
    )

  # =====> CONVERT TO RANDOM INTERCEPT?
  # Where prior_def is for sd, where randint ~ N(0, sd)
  if (param == "sigma") {
    prior_samples <- rnorm(nsamples, mean = 0, sd = prior_samples)
    plots$random_intercept <- .plot_samples(prior_samples, ndots, point_interval, interval_widths, interval_colors, na_color) + 
      labs(
        x = latex2exp::TeX("\\Delta"), 
        title =  latex2exp::TeX("\\Delta ~ Normal(0, \\sigma)")
      )
    message("sd prior converted to random intercept as RI ~ N(0, sd)")
  }
  # =====> BACKTRANSFORM TO RESPONSE SCALE?
  if (link != "identity") {
    inv_link <- .inverse_link(link)
    prior_samples <- .backtransform(prior_samples, inv_link)
    
    if (param == "sigma") {
      xlab <- paste0(inv_link, "(\\Delta)")
      vline_col <- vline_color
    } else {
      vline_col <- NULL
      xlab <- paste0(inv_link, "(\\", param, ")")
    }
    plots$backtransformed <- .plot_samples(
      prior_samples, ndots, point_interval, 
      interval_widths, interval_colors, na_color, vline_color = vline_col
    ) + 
    labs(x = latex2exp::TeX(xlab)) 
    message(paste0("inverse link function applied [link: ", link, "; inverse link: ", inv_link), "]")
  }
  # =====> RETURN PLOTS
  if (patched) {
    plots |> 
      patchwork::wrap_plots(ncol = 1) + 
      patchwork::plot_layout(guides = 'collect')
  } else {
    return(plots)
  }
}


# ==================================
#         HELPER FUNCTIONS
# ==================================
.plot_samples <- function(samples, ndots, point_interval, interval_widths, interval_colors, na_color, vline_color = NULL) {

  # NOTE: ggdist seems to be unreliable in its point calculations, e.g. 
  # set.seed(710); mode_hdci(exp(rnorm(1e5, 0, rexp(1e5, 2)))
  # 
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
  
  # CREATE PLOT
  p <- data.frame(x = samples) |> ggplot()
  if (!is.null(vline_color)) p <- p + geom_vline(xintercept = 1, linetype = "dotted", color = vline_color)
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
      slab_fill = "HDCI", slab_color = "HDCI"
    ) + 
    gg_remove_yaxis

  return(p)
}
.get_prior_samples <- function(prior_def, n_samples, seed, msg = FALSE) {
  synonyms <- list(
    rnorm = c("normal", "gaussian", "N"), 
    rgamma = c("gamma"),
    rexp = c("exponential", "exp"), 
    runif = c("uniform", "U", "unif")
  )
  prior_def <- str_split_1(prior_def, "~")[2] |> str_trim() # e.g. beta ~ N(0, 1) -> N(0, 1)
  distr <- prior_def |> 
    str_extract("^\\w+") |>  # e.g. N(0, 1) -> N
    str_to_lower()
  rfunc <- names(which(sapply(synonyms, function(x) distr %in% tolower(x)))) # e.g. N -> rnorm
  if (nchar(rfunc) == 0) stop(paste0("prior_def must be one of [", paste(unlist(unlist(synonyms)), collapse = ", "), "]. Add distributions via `synonyms` in get_prior_samples()"))
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
  if (!(length(prior_samples) == n_samples)) stop(paste0("Trying to call `", rfunc_call, "` does not generate the expected number of samples."))
  if (msg) message(paste0(length(prior_samples), " prior samples produced via `", rfunc_call, "`"))
  return(prior_samples)
}
.map_samples_to_hdci_intervals <- function(prior_samples, hdci) {
  hdci <- hdci |> arrange(desc(.width))
  which_hdci <- rep(NA, length(prior_samples))
  for (i in 1:nrow(hdci)) {
    width_i <- as.numeric(hdci[i, ".width"])
    lower_i <- as.numeric(hdci[i, "ymin"])
    upper_i <- as.numeric(hdci[i, "ymax"])
    which_hdci[prior_samples >= lower_i & prior_samples <= upper_i] <- width_i
  }
  return(which_hdci)
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
.get_param_name <- function(prior_def) {
  if (!str_detect(prior_def, "~")) stop("prior_def should be of form '<parameter> ~ <distr_name>(<distr_parameters>)'")
  param <- str_trim(str_split_1(prior_def, "~")[1])
  p <- tolower(param)
  if (p %in% c("sd", "sigma")) param <- "sigma"
  if (param %in% c("b", "beta")) param <- "beta"
  return(param)
}