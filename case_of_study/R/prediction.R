#'  Prediction_BART function
#'  
#' Make Predictions Using BART Model This function makes predictions using a BART model
#' on a stack of environmental covariates.
#' @param bart_model A BART model.
#' @param stack_cov A raster stack containing environmental covariates for prediction.
#' @return A SpatRaster containing the mean, median, standard deviation, and quantiles of the posterior predictive distribution 
#' @export
prediction_BART <- function(bart_model, stack_cov) {
  tryCatch({
    quantiles <- c(0.025, 0.975)
    input.matrix <- terra::as.matrix(stack_cov)
    blankout <- data.frame(matrix(ncol=(4+length(quantiles)), 
                                  nrow=ncell(stack_cov[[1]])))
    whichvals <- which(complete.cases(input.matrix))
    input.matrix <- input.matrix[complete.cases(input.matrix),]
    
    pred <- dbarts:::predict.bart(bart_model, input.matrix)
    
    pred_data <- cbind(data.frame(matrixStats::colMeans2(pred)),
                       data.frame(matrixStats::colMedians(pred)),
                       data.frame(matrixStats::colSds(pred)),
                       data.frame(matrixStats::colQuantiles(pred, probs = quantiles)))
    
    pred_data$diff <- pred_data$X97.5.-pred_data$X2.5.
    names(pred_data) <- c("mean","median","sd","q0.025","q0.975","diff")
    
    blankout[whichvals,] <- as.matrix(pred_data)
    output <- blankout
    
    outlist <- lapply(1:ncol(output), function(x) {
      output.m <- t(matrix(output[,x],
                           nrow = ncol(stack_cov),
                           ncol = nrow(stack_cov)))
      return(terra::rast(output.m, extent = c(-180, 180, -90, 90)))
    })
    
    outlist <- terra::rast(outlist)
    names(outlist) <- c("mean", "median", "sd", "q0.025", "q0.975", "diffquantiles")
    
    return(outlist)
  }, error = function(err) {
    message("Model failed: ", conditionMessage(err))
  })
}
