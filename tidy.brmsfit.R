tidy.brmsfit <- function(x, ...) {
  sfixed <- summary(x)$fixed
  fixed <- data.frame()
  if (!is.null(sfixed)) {
    fixed <- data.frame(
      effect = rep("fixed", nrow(sfixed)),
      group = rep("NA", nrow(sfixed)),
      term = rownames(sfixed), 
      estimate = sfixed$Estimate,
      std.error = sfixed$Est.Error,
      conf.low = sfixed$`l-95% CI`,
      conf.high = sfixed$`u-95% CI`
    )
  }
  srandom <- summary(x)$random
  random <- data.frame()
  if (!is.null(srandom)) {
    for (group in names(srandom)) {
      grandom <- srandom[[group]]
      random_group <- data.frame(
        effect = rep("random", nrow(grandom)),
        group = rep(group, nrow(grandom)),
        term = rownames(grandom),
        estimate = grandom$Estimate,
        std.error = grandom$Est.Error, 
        conf.low = grandom$`l-95% CI`,
        conf.high = grandom$`u-95% CI`
      )
      random <- rbind(random, random_group)
    } 
  }
  rbind(fixed, random)
}