#' Optimal cutoff for presence-absence prediction
#'
#' This function calculates the optimal cutoff for presence-absence prediction using a BART model.
#'
#' @param pa_coords Dataframe with coordinates of presences and absences and column indicating if presence or absence.
#' @param layers A list of layer names to extract from the raster stack.
#' @param model A BART model object.
#' @param seed Random seed for reproducibility.
#'
#' @return The optimal cutoff value for presence-absence prediction.
#'
#' @export
pa_optimal_cutoff <- function(pa_coords, layers, model, seed = NULL) {
  data <- terra::extract(layers, pa_coords[, c("decimalLongitude", "decimalLatitude")])
  set.seed(seed)
  pred_data <- dbarts:::predict.bart(model, newdata = data[, names(layers), drop = FALSE])
  pred_mean <- colMeans(pred_data)
  pa_cutoff <- optimalCutoff(
    actuals = pa_coords[, "pa"],
    predictedScores = pred_mean
  )
  return(pa_cutoff)
}
#' Extract Covariate Values without NA
#'
#' This function extracts covariate values for species occurrences, excluding NA values.
#'
#' @param data A data frame containing species occurrence data with columns "decimalLongitude" and "decimalLatitude".
#' @param covariate_layers A list of raster layers representing covariates.
#' @return A data frame containing species occurrence data with covariate values, excluding NA values.
#' @details This function extracts covariate values for each species occurrence location from the provided covariate layers. It returns a data frame containing species occurrence data with covariate values, excluding any NA values.
#'
#' @export
extract_noNA_cov_values <- function(data, covariate_layers){
  covariate_values <- terra::extract(
    covariate_layers,
    data[, c("decimalLongitude", "decimalLatitude")]
  )
  covariate_values <- cbind(data, covariate_values) %>%
    tidyr::drop_na()
  return(covariate_values)
}

cv_bart <- function(pa_coords, layers, k = 5, seed = NULL){
  set.seed(seed)
  # Extract covariate values
  pa_coords <- extract_noNA_cov_values(pa_coords, layers)
  n <- nrow(pa_coords)
  # Create index vector
  k_index <- rep(1:k, length.out = n)
  # Randomize vector
  k_index <- sample(k_index)
  pa_coords$k <- k_index
  TP <- c()
  FP <- c()
  FN <- c()
  TN <- c()
  for (i in 1:k){
    # Split train-test
    train <- pa_coords[pa_coords$k != i, ]
    test <- pa_coords[pa_coords$k == i, ]
    bart_model <- dbarts::bart(x.train = train[, names(layers), drop = FALSE],
                               y.train = train[,"pa"],
                               keeptrees = TRUE)
    invisible(bart_model$fit$state)
    cutoff <- pa_optimal_cutoff(train, layers, bart_model)
    pred <- dbarts:::predict.bart(bart_model, test[, names(layers), drop = FALSE])
    pred <- colMeans(pred)
    potential_presences <- ifelse(pred >= cutoff, 1, 0)
    # Confusion matrices
    TP <- c(TP, sum(test$pa == 1 & potential_presences == 1))
    FP <- c(FP, sum(test$pa == 0 & potential_presences == 1))
    FN <- c(FN, sum(test$pa == 1 & potential_presences == 0))
    TN <- c(TN, sum(test$pa == 0 & potential_presences == 0))
  }
  # Metrics
  PREC = TP / (TP + FP)
  SEN = TP / (TP + FN)
  SPC = TN / (TN + FP)
  FDR = FP / (TP + FP)
  NPV = TN / (FN + TN)
  FNR = FN / (TP + FN)
  FPR = FP / (FP + TN)
  Fscore = 2 * ((PREC * SEN) / (PREC + SEN))
  ACC = (TP + TN)/(TP + FP + FN + TN)
  BA = (SEN + SPC) / 2
  TSS = SEN + SPC - 1
  cv_res <- data.frame(
    TP = TP, FP = FP, FN = FN, TN = TN, PREC = PREC, SEN = SEN, SPC = SPC,
    Fscore = Fscore, ACC = ACC, BA = BA, FNR = FNR, FPR = FPR, FDR = FDR,
    NPV = NPV, TSS = TSS
  )
  return(cv_res)
}
