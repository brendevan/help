glance.brmsfit <- function(x, ...) {
  s <- summary(x)
  r2 <- performance::r2_bayes(x) |> as.data.frame()
  if (nrow(r2) == 1) {
    r2.cond <- r2$R2
    r2.cond.sd <- r2$SD
    r2.cond.95CI.low <- r2$CI_low
    r2.cond.95CI.high <- r2$CI_high
    r2.marg <- NA_real_
    r2.marg.sd <- NA_real_
    r2.marg.95CI.low <- NA_real_
    r2.marg.95CI.high <- NA_real_
  } else {
    r2c <- dplyr::filter(r2, Component == "conditional")
    r2.cond <- r2c$R2
    r2.cond.sd <- r2c$SD
    r2.cond.95CI.low <- r2c$CI_low
    r2.cond.95CI.high <- r2c$CI_high
    r2m <- dplyr::filter(r2, Component == "marginal")
    r2.marg <- r2m$R2
    r2.marg.sd <- r2m$SD
    r2.marg.95CI.low <- r2m$CI_low
    r2.marg.95CI.high <- r2m$CI_high
  }
  loo <- loo(x)
  pareto_k <-  loo$diagnostics$pareto_k
  loo_est <- loo$estimates |> as.data.frame()
  waic <- brms::waic(x)$estimates
  interval_count <- function(x, from, to) sum(x > from & x <= to)
  data.frame(
    n = s$nobs,
    # R squared 
    r2.cond, r2.cond.sd, r2.cond.95CI.low, r2.cond.95CI.high,
    r2.marg, r2.marg.sd, r2.marg.95CI.low, r2.marg.95CI.high, 
    # WAIC
    elpd.waic = waic["elpd_waic", "Estimate"], 
    se.elpd.waic = waic["elpd_waic", "SE"], 
    p.waic = waic["p_waic", "Estimate"], 
    se.p.waic = waic["p_waic", "SE"], 
    waic = waic["waic", "Estimate"],
    se.waic = waic["elpd_waic", "SE"], 
    # Leave-one-out metrics
    elpd.loo = loo_est["elpd_loo", "Estimate"], 
    se.elpd.loo = loo_est["elpd_loo", "SE"], 
    p.loo = loo_est["p_loo", "Estimate"], 
    se.p.loo = loo_est["p_loo", "SE"], 
    looic = loo_est["looic", "Estimate"], 
    se.looic = loo_est["looic", "SE"], 
    # Leave-one-out diagnostics (also apply to WAIC)
    pareto.k.ltp5 = interval_count(pareto_k, -Inf, 0.5), 
    pareto.k.p5top7 = interval_count(pareto_k, 0.5, 0.7), 
    pareto.k.p7to1 = interval_count(pareto_k, 0.7, 1), 
    pareto.k.gt1 = interval_count(pareto_k, 1, Inf)
  )
}