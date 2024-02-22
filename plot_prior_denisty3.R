
plot_prior_density3 <- function(
  prior_def, # Prior distribution definition used to generate samples
  n_samples = 1000, # Number of samples to generate
  hdci_widths = c(0.5, 0.89, 0.99), # HDCI widths to plot
  hdci_colors = c("#FF6348FF", "#008EA0FF", "#8A4198FF"),
  subtitle_width = 60, # Number of characters after which the subtitle is wrapped to the next line
  adjust_margin = FALSE, # Try to fix the top margin of the plot if subtitle is not fully shown, 
  hdci_table = "none",    # Whether or not to print the table of HDCIs: "none", "print", "return"\
  seed = "random",       # Integer; if "random" (default) a random integer in [0, 1000] is used
  ...                    # Arguments for ggdist::geom_dotsinterval(); e.g. alpha = 0.5
) {
  require(ggplot2)
  require(ggdist)
  require(ggtext)
  require(dplyr)
  require(stringr)

  # =====> CHECKS
  if (length(hdci_widths) != length(hdci_colors)) {
    stop("hdci_widths and hdci_colors must be the same length")
  }
  # =====> GET PRIOR SAMPLES AND HDCI
  if (seed == "random") seed <- round(runif(1)*1000, 0)
  message(paste("Random seed =", seed))
  prior_samples <- .get_prior_samples(prior_def, n_samples, seed)
  message(paste0(length(prior_samples), " prior samples produced via `", rfunc_call, "`"))
  million_prior_samples <- .get_prior_samples(prior_def, 1e6, seed)
  hdci <- million_prior_samples |> 
    ggdist::mode_hdci(.width = hdci_widths)
  if (hdci_table == "print") print(hdci)
  if (hdci_table == "return") return(hdci)
  prior_samples <- data.frame(prior_sample = prior_samples) |> 
    mutate(
      hdci = .map_samples_to_hdci_intervals(prior_sample, hdci), 
      hdci_labels = case_when(
        is.na(hdci) ~ "", 
        TRUE ~ paste0(round(hdci*100, 0), "%")  
      )
    )
  prior_samples |> 
    ggplot() + 
    geom_dots(aes(x = prior_sample, fill = hdci_labels), slab_color = NA, alpha = 0.8) + 
    # geom_interval(data = hdci, aes(y = 0, x = y, xmin = ymin, xmax = ymax, thickness = c(1, 2, 3))) + 
    scale_fill_manual(values = hdci_colors, na.value = "black") + 
    labs(fill = "HDCI", x = "Prior samples") 
}

.get_prior_samples <- function(prior_def, n_samples, seed) {
  synonyms <- list(
    rnorm = c("normal", "gaussian", "N"), 
    rgamma = c("gamma"),
    rexp = c("exponential", "exp"), 
    runif = c("uniform", "U", "unif")
  )
  distr <- tolower(str_extract(prior_def, "^\\w+"))
  rfunc <- names(which(sapply(synonyms, function(x) distr %in% tolower(x))))
  if (nchar(rfunc) == 0) stop(paste0("prior_def must be one of [", paste(unlist(unlist(synonyms)), collapse = ", "), "]. Add distributions via `synonyms` in get_prior_samples()"))
  args <- tolower(prior_def) |> 
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
  prior_samples
}
.map_samples_to_hdci_intervals <- function(prior_samples, hdci) {
  hdci <- hdci |> arrange(desc(hdci_width))
  which_hdci <- rep(NA, length(prior_samples))
  for (i in 1:nrow(hdci)) {
    width_i <- as.numeric(hdci[i, "hdci_width"])
    lower_i <- as.numeric(hdci[i, "lower"])
    upper_i <- as.numeric(hdci[i, "upper"])
    which_hdci[prior_samples >= lower_i & prior_samples <= upper_i] <- width_i
  }
  return(which_hdci)
}
