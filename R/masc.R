#' MASC - Mixed effect modeling of Associations of Single Cells
#'
#' @param dataset A data frame containing the contrast factor, random, and fixed effects for the model
#' @param cluster A factor indicating cluster assignments for each cell
#' @param contrast A vector indicating the variable to be tested for association with cluster abundance. Must match a column in dataset.
#' @param random_effects A vector indicating which terms should be modeled as random effects covariates. Terms listed must match columns in dataset.
#' @param fixed_effects A vector indicating which terms should be modeled as fixed effects covariates. Terms listed must match columns in dataset.
#' @param save_models Should MASC save the mixed-effects model objects generated for each cluster?
#' @param save_model_dir Location to save mixed-effect model objects. Defaults to current working directory.
#' @param verbose TRUE/FALSE
#'
#' @return data frame containing calculated association p-values and odds ratios for each cluster tested
#'
#' @examples
#' # Create test dataset with three clusters of 100 cells each
#' test.df <- data.frame(cluster = factor(rep(c(1, 2, 3), each = 100)))
#' # Create 6 donors that are cases or controls and include covariates
#' donors.df <- data.frame(donor = rep(paste("Donor", LETTERS[1:6], sep = "_"), each = 50),
#' sex = rep(c("M", "F", "M", "F", "F", "M"), each = 50),
#' status = rep(c("Case", "Case", "Control", "Control", "Case", "Control"), each = 50))
#' # Now make cluster 1 mostly case, cluster 2 mostly controls, etc
#' cases <- donors.df[donors.df$status == "Case",]
#' cases <- cases[sample(nrow(cases)),]
#' controls <- donors.df[donors.df$status == "Control",]
#' controls <- controls[sample(nrow(controls)),]
#' test.df <- cbind(rbind(cases[1:75,], controls[1:25,], cases[76:115,], controls[26:85,], cases[116:150,], controls[86:150,]), test.df)
#' # Test set call
#' library(lme4)
#' MASC(data = test.df, cluster = test.df$cluster, contrast = "status", random.effects = "donor", fixed.effects = "sex")
#'

MASC <- function(dataset, cluster, contrast, random_effects = NULL, fixed_effects = NULL,
                 verbose = FALSE, save_models = FALSE, save_model_dir = NULL) {
  # Check inputs
  if (is.factor(dataset[[contrast]]) == FALSE) {
    stop("Specified contrast term is not coded as a factor in dataset")
  }

  # Generate design matrix from cluster assignments
  cluster <- as.character(cluster)
  designmat <- model.matrix(~ cluster + 0, data.frame(cluster = cluster))
  dataset <- cbind(designmat, dataset)

  # Convert cluster assignments to string
  cluster <- as.character(cluster)
  # Prepend design matrix generated from cluster assignments
  designmat <- model.matrix(~ cluster + 0, data.frame(cluster = cluster))
  dataset <- cbind(designmat, dataset)
  # Create output list to hold results
  res <- vector(mode = "list", length = length(unique(cluster)))
  names(res) <- attributes(designmat)$dimnames[[2]]

  # Create model formulas
  if (!is.null(fixed_effects) && !is.null(random_effects)) {
    model_rhs <- paste0(c(paste0(fixed_effects, collapse = " + "),
                          paste0("(1|", random_effects, ")", collapse = " + ")),
                        collapse = " + ")
    if (verbose == TRUE) {
      message(paste("Using null model:", "cluster ~", model_rhs))
    }
  } else if (!is.null(fixed_effects) && is.null(random_effects)) {
    model_rhs <- paste0(fixed_effects, collapse = " + ")
    if (verbose == TRUE) {
      message(paste("Using null model:", "cluster ~", model_rhs))
      # For now, do not allow models without mixed effects terms
      stop("No random effects specified")
    }
  } else if (is.null(fixed_effects) && !is.null(random_effects)) {
    model_rhs <- paste0("(1|", random_effects, ")", collapse = " + ")
    if (verbose == TRUE) {
      message(paste("Using null model:", "cluster ~", model_rhs))
    }
  } else {
    model_rhs <- "1" # only includes intercept
    if (verbose == TRUE) {
      message(paste("Using null model:", "cluster ~", model_rhs))
      stop("No random or fixed effects specified")
    }
  }


  # Initialize list to store model objects for each cluster
  cluster_models <- vector(mode = "list",
                           length = length(attributes(designmat)$dimnames[[2]]))
  names(cluster_models) <- attributes(designmat)$dimnames[[2]]

  # Run nested mixed-effects models for each cluster
  for (i in seq_along(attributes(designmat)$dimnames[[2]])) {
    test_cluster <- attributes(designmat)$dimnames[[2]][i]
    if (verbose == TRUE) {
      message(paste("Creating logistic mixed models for", test_cluster))
    }
    null_fm <- as.formula(paste0(c(paste0(test_cluster, " ~ 1 + "),
                                   model_rhs), collapse = ""))
    full_fm <- as.formula(paste0(c(paste0(test_cluster, " ~ ", contrast, " + "),
                                   model_rhs), collapse = ""))
    # Run null and full mixed-effects models
    null_model <- lme4::glmer(formula = null_fm, data = dataset,
                              family = binomial, nAGQ = 1, verbose = 0,
                              control = glmerControl(optimizer = "bobyqa"))
    full_model <- lme4::glmer(formula = full_fm, data = dataset,
                              family = binomial, nAGQ = 1, verbose = 0,
                              control = glmerControl(optimizer = "bobyqa"))
    model_lrt <- anova(null_model, full_model)
    # calculate confidence intervals for contrast term beta
    contrast_lvl2 <- paste0(contrast, levels(dataset[[contrast]])[2])
    contrast_ci <- confint.merMod(full_model, method = "Wald",
                                  parm = contrast_lvl2)
    # Save model objects to list
    cluster_models[[i]]$null_model <- null_model
    cluster_models[[i]]$full_model <- full_model
    cluster_models[[i]]$model_lrt <- model_lrt
    cluster_models[[i]]$confint <- contrast_ci
  }

  # Organize results into output dataframe
  output <- data.frame(cluster = attributes(designmat)$dimnames[[2]],
                       size = colSums(designmat))
  output$model.pvalue <- sapply(cluster_models, function(x) x$model_lrt[["Pr(>Chisq)"]][2])
  output[[paste(contrast_lvl2, "OR", sep = ".")]] <- sapply(cluster_models, function(x) exp(fixef(x$full)[[contrast_lvl2]]))
  output[[paste(contrast_lvl2, "OR", "95pct.ci.lower", sep = ".")]] <- sapply(cluster_models, function(x) exp(x$confint[contrast_lvl2, "2.5 %"]))
  output[[paste(contrast_lvl2, "OR", "95pct.ci.upper", sep = ".")]] <- sapply(cluster_models, function(x) exp(x$confint[contrast_lvl2, "97.5 %"]))

  # Return MASC results and save models if specified
  if (save_models == TRUE) {
    saveModelObj(cluster_models, save_dir = save_model_dir)
    return(output)
  } else {
    return(output)
  }
}



