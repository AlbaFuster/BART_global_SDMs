cov_relations <- function(bart_model, data, names) {
  
  values_q <- list()
  
  for (i in 1:length(names)) {
    column_values <- na.omit(data[, names[i]])
    quantiles <- quantile(column_values, prob = c(0.05, seq(0.1, 0.9, 0.05), 0.95))
    
    if (any(duplicated(quantiles))) {
      warning(paste("Duplicate values found in quantiles for column", names[i]))
    }
    
    values_q[[i]] <- quantiles
  }
  
  var_r <- dbarts::pdbart(
    bart_model,
    xind = names,
    levsquants = c(0.05, seq(0.1, 0.9, 0.05), 0.95)
  )
  
  inv_probit_df <- list()
  for(i in 1:length(names)) {
    inv_probit_df[[i]] <- data.frame(matrix(nrow = nrow(var_r$fd[[i]]), ncol = dim(var_r$fd[[i]])[2])) 
  }
  
  for(j in 1:length(names)) {
    for (i in 1:dim(var_r$fd[[1]])[2]) {
      inv_probit_df[[j]][, i] <- mcp::phi(var_r$fd[[j]][,i])
    } 
  }
  
  prob_inv <- list()
  quantile_values_q25 <- list()
  quantile_values_q975 <- list()
  values_var <- list()
  
  # Loop through the elements in the 'names' list
  for (i in seq_along(names)) {
    # Calculate mean and store in prob_inv
    prob_inv[[i]] <- as.vector(rbind(sapply(inv_probit_df[[i]][, 1:dim(var_r$fd[[i]])[2]], function(x) mean(x))))
    
    # Calculate the 2.5th percentile and store in quantile_values_q25
    quantile_values_q25[[i]] <- as.vector(apply(inv_probit_df[[i]][, 1:dim(var_r$fd[[i]])[2]], 2, function(x) quantile(x, probs = 0.025, na.rm = TRUE)))
    
    # Calculate the 97.5th percentile and store in quantile_values_q975
    quantile_values_q975[[i]] <- as.vector(apply(inv_probit_df[[i]][, 1:dim(var_r$fd[[i]])[2]], 2, function(x) quantile(x, probs = 0.975, na.rm = TRUE)))
    
    # Store values from var_r$levs in values_var
    values_var[[i]] <- var_r$levs[[i]]
  }
  
  
  data_var <- list()
  for(i in 1:length(names)) {
    data_var[[i]] <- data.frame(prob = prob_inv[[i]],
                                q25 = quantile_values_q25[[i]],
                                q975 = quantile_values_q975[[i]],
                                var = values_var[[i]]
    ) 
  }
  
  return(data_var)
}
