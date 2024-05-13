#' @title Variable importance plot
#'
#' @description
#'
#' Variable importance, as measured in the proportion of total branches used for a given variable.
#' 
#' @param bart_model A BART model. 
#' @return A numeric vector containing the name of the variable and its corresponding importance value. The variables are sorted such that the first variable listed is the one that has the greatest influence.
#'
varimp <- function(bart_model) {
  importance <- colMeans(bart_model$varcount/rowSums(bart_model$varcount))
  importance <- sort(importance, decreasing = TRUE)
  return(importance)
}