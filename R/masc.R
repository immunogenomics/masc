#' MASC - Mixed effect modeling of Associations of Single Cells
#'
#' @param data A data frame containing the contrast factor, random, and fixed effects for the model
#' @param cluster A factor indicating cluster assignments for each cell
#' @param contrast A vector indicating the factor to be used as a contrast term
#' @param random.effects A vector indicating which terms should be modeled as random effects
#' @param fixed.effects A vector indicating wich terms should be modeled as fixed effects
#' @param verbose TRUE/FALSE
#'

MASC <- function(data, cluster, contrast, random.effects, fixed.effects, verbose = TRUE) {
  # Create design matrix from cluster assignments
  M <- model.matrix(~ cluster + 0, data.frame(cluster = cluster))
  # Combine into data frame with other effect terms
  df <- cbind(as.data.frame(M), data)
  # Create output list to hold results
  modelList <- vector(mode = "list", length = nlevels(cluster))
  names(modelList) <- colnames(M)
  # Run mixed effect model on each cluster

  # Create model formulas
  model.randeff <- paste("(1|", random.effects, ")", sep = "", collapse = " + ")
  model.fixeff <- paste(fixed.effects, collapse = " + ")
  for (i in 1:nlevels(cluster)) {
    clusterName <- levels(cluster)[i]
    message(paste("Creating logistic mixed models for", clusterName))
    null.model <- as.formula(paste(paste("cluster", clusterName, " ~ 1", sep = ""), model.fixeff, model.randeff, sep = " + "))
    full.model <- as.formula(paste(paste("cluster", clusterName, " ~ ", contrast, sep = ""), model.fixeff, model.randeff, sep = " + "))
    modelList[[i]]$null <- glmer(formula = null.model, data = df, family = binomial,
                                 control = glmerControl(optimizer = "bobyqa") ,nAGQ = 1,  verbose = 0)
    modelList[[i]]$full <- glmer(formula = full.model, data = df, family = binomial,
                                 control = glmerControl(optimizer = "bobyqa") ,nAGQ = 1,  verbose = 0)
    modelList[[i]]$anova <- anova(modelList[[i]]$null, modelList[[i]]$full)
    # calculate confidence intervals for the case-control beta
    contrast.ci <- paste(contrast, levels(data[[contrast]])[2], sep = "")
    modelList[[i]]$confint <- confint.merMod(modelList[[i]]$full, parm = contrast.ci, method = "Wald", devtol = 1e-6)
  }
  out.df <- data.frame(cluster = colnames(M), cells = colSums(M))
  # Add logistic mixed-effect model p-values and odds ratios for each cluster
  out.df$masc.pval <- sapply(modelList, function(x) x$anova[["Pr(>Chisq)"]][2])
  out.df$masc.or <- sapply(modelList, function(x) exp(fixef(x$full)[[contrast.ci]]))
  out.df$masc.or.lower <- sapply(modelList, function(x) exp(x$confint[contrast.ci, "2.5 %"]))
  out.df$masc.or.upper <- sapply(modelList, function(x) exp(x$confint[contrast.ci, "97.5 %"]))
  return(out.df)
}


