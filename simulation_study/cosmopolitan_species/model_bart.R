#*********************************
# BART with simulated data       #
# Biomass and sampling           #
#*********************************
# Alba Fuster-Alonso             #
#********************************#

library(dbarts)
library(embarcadero)
library(raster)
library(caret)
library(Metrics)

load("./simulated_data_sampling.RData")

# Response variable
data_random_17 <- list()
for(i in 1:50) {
  data_random_17[[i]] <- data_random_final[[i]][!(data_random_final[[i]]$time >= 18 & data_random_final[[i]]$time <= 20), ]
}

#########################
#FITTING#################
#########################
bart_fit <- list()
for(i in 1:50) {
  bart_fit[[i]] <- dbarts::bart(x.train = data_random_17[[i]][,c(2,3,4,5,7)],
                           y.train = data_random_17[[i]]$ocurrences,
                           keeptrees = TRUE)
}

for(i in 1:50) {
  print(i)
  invisible(bart_fit[[i]]$fit$state) # very important if you want to use the models later
}

for(i in 1:50) {
  summary(bart_fit[[i]]) # extract cutoff!! NO LO HE CONSEGUIDO
}

cutoff <- c(0.4865212,0.4953986,0.6310514,0.5910016,0.5207975,0.6067206,0.5966403,
            0.611047,0.4416954,0.4644592,0.4868417,0.5259708,0.51204,0.4980275,
            0.5664022,0.5148826,0.560783,0.4517209,0.5546841,0.5240254,0.5697541,
            0.555204,0.4512806,0.4872573,0.6101182,0.4909938,0.5491327,0.5234882,
            0.5835319,0.5563848,0.4502762,0.4478078,0.4855068,0.5546058,0.5734038,
            0.4474298,0.4986877,0.4424279,0.4179429,0.5916037,0.4766743,0.5820038,
            0.462559,0.5072233,0.5647027,0.5703927,0.4736139 ,0.4673741,0.5072344,
            0.5631959)

fitted_values <- list()
for(i in 1:50) {
  fitted_values[[i]] <- fitted(bart_fit[[i]])
  data_random_17[[i]]$fit <- fitted_values[[i]]
  data_random_17[[i]]$fit_ocurrence <- ifelse(data_random_17[[i]]$fit >= cutoff[i], 1, 0)
}

# Calculate measures
conf_matrix <- list()
for(i in 1:50) {
  conf_matrix[[i]] <- table(data_random_17[[i]]$ocurrences,
                            data_random_17[[i]]$fit_ocurrence) 
}

sensitivity_intra <- list()
specificity_intra <- list()
accuracy_intra <- list()
for(i in 1:50) {
  sensitivity_intra[[i]] <- sensitivity(conf_matrix[[i]])
  specificity_intra[[i]] <- specificity(conf_matrix[[i]])
  accuracy_intra[[i]] <- accuracy(data_random_17[[i]]$ocurrences,
                                  data_random_17[[i]]$fit_ocurrence)
}

data_measures_intra <- data.frame(
  sensitivity = unlist(sensitivity_intra),
  specificity = unlist(specificity_intra),
  accuracy = unlist(accuracy_intra)
)

#########################
#PREDICTION##############
#########################

# Variables
for(i in 1:k) {
  data_variables[[i]]$time <- rep(i, 10000) 
}

sim_bat <- rasterFromXYZ(data_variables[[1]][,c(4,5,2)]) 
sim_xcoord <- rasterFromXYZ(data_variables[[1]][,c(4,5,4)])
sim_ycoord <- rasterFromXYZ(data_variables[[1]][,c(4,5,5)])

sim_year <- list()
for(i in 1:k) {
  sim_year[[i]] <- rasterFromXYZ(data_variables[[i]][,c(4,5,7)]) 
}

sim_temp <- list()
for(i in 1:k) {
  sim_temp[[i]] <- rasterFromXYZ(data_variables[[i]][,c(4,5,3)])
}

for(i in 1:k) {
  writeRaster(sim_bat,
              filename = paste0("./covariate/bat",2000 + i,".tif"),
              options=c('TFW=NO'))
}

for(i in 1:k) {
  writeRaster(sim_xcoord,
              filename = paste0("./covariate/xcoord",2000 + i,".tif"),
              options=c('TFW=NO'))
}

for(i in 1:k) {
  writeRaster(sim_ycoord,
              filename = paste0("./covariate/ycoord",2000 + i,".tif"),
              options=c('TFW=NO'))
}

for(i in 1:k) {
  writeRaster(sim_year[[i]],
              filename = paste0("./covariate/year",2000 + i,".tif"),
              options=c('TFW=NO'))
}

for(i in 1:k) {
  writeRaster(sim_temp[[i]],
              filename = paste0("./covariate/temp",2000 + i,".tif"),
              options=c('TFW=NO'))
}


years <- 2001:2020
raster_files <- list()
for(i in 1:length(years)) {
  raster_files[[i]] <- list.files(path = "./covariate/",
                                                pattern = paste0(as.character(years[i])))
}

raster_list <- list()
ref_raster <- NULL  # Raster de referencia para asegurarse de que todos los rasters tengan la misma extensión

for (i in seq_along(raster_files)) {
  print(i)
  stack_list <- list()
  for (j in seq_along(raster_files[[i]])) {
    print(j)
    file_path <- paste0("./covariate/", raster_files[[i]][j])
    raster <- raster::raster(file_path)
    if (is.null(ref_raster)) {  # Si es el primer raster, usarlo como referencia
      ref_raster <- raster
    } else if (!compareRaster(raster, ref_raster, stopiffalse = FALSE)) {  # Si la extensión es diferente, hacer una reproyección
      raster <- projectRaster(raster, ref_raster)
    }
    stack_list[[j]] <- raster
  }
  raster_stack <- do.call(raster::stack, stack_list)
  raster_list[[i]] <- raster_stack
}

for(i in 1:k) {
  names(raster_list[[i]]) <- c("batimetria","temp","xcoord","ycoord","time")
}

pred_BART <- list()
for (j in 1:50) {
  print(j)
  pred_BART_j <- list()  # Nueva lista para cada iteración de j
  for (i in 1:k) {
    print(i)
    pred_BART_j[[i]] <- predict2.bart(bart_fit[[j]],
                                      x.layers = raster_list[[i]],
                                      quantiles = c(0.025, 0.5, 0.975)) 
    names(pred_BART_j[[i]]) <- c("predicted", "q0025", "median", "q0975")
  }
  pred_BART[[j]] <- pred_BART_j
}

save.image(file = "./fit_pred.RData")

prediction_all <- list()
for(j in 1:50) {
  print(j)
  pred_BART_ext <- list()
  for (i in 1:20) {
    print(i)
    coordinates <- raster::xyFromCell(pred_BART[[j]][[i]]$median, 1:ncell(pred_BART[[j]][[i]]$median))
    median_values <- raster::extract(pred_BART[[j]][[i]]$median, coordinates)
    
    df <- data.frame(coordinates, median = median_values)
    pred_BART_ext[[i]] <- df
  }
  prediction_all[[j]] <- pred_BART_ext
}

prediction_order <- list()
for(j in 1:50) {
  prediction_order_1 <- list()
  for (i in 1:k) {
    prediction_order_1[[i]] <- prediction_all[[j]][[i]][order(prediction_all[[j]][[i]]$y),]
  } 
  prediction_order[[j]] <- prediction_order_1
}
                                                                                                           

data_all_2 <- list()
for(j in 1:50) {
  print(j)
  data_all_2[[j]] <- data.frame(probability = c(data_variables[[1]]$probability,
                                         data_variables[[2]]$probability,
                                         data_variables[[3]]$probability,
                                         data_variables[[4]]$probability,
                                         data_variables[[5]]$probability,
                                         data_variables[[6]]$probability,
                                         data_variables[[7]]$probability,
                                         data_variables[[8]]$probability,
                                         data_variables[[9]]$probability,
                                         data_variables[[10]]$probability,
                                         data_variables[[11]]$probability,
                                         data_variables[[12]]$probability,
                                         data_variables[[13]]$probability,
                                         data_variables[[14]]$probability,
                                         data_variables[[15]]$probability,
                                         data_variables[[16]]$probability,
                                         data_variables[[17]]$probability,
                                         data_variables[[18]]$probability,
                                         data_variables[[19]]$probability,
                                         data_variables[[20]]$probability),
                         
                         pred_probability = c(
                           prediction_order[[j]][[1]]$median,
                           prediction_order[[j]][[2]]$median,
                           prediction_order[[j]][[3]]$median,
                           prediction_order[[j]][[4]]$median,
                           prediction_order[[j]][[5]]$median,
                           prediction_order[[j]][[6]]$median,
                           prediction_order[[j]][[7]]$median,
                           prediction_order[[j]][[8]]$median,
                           prediction_order[[j]][[9]]$median,
                           prediction_order[[j]][[10]]$median,
                           prediction_order[[j]][[11]]$median,
                           prediction_order[[j]][[12]]$median,
                           prediction_order[[j]][[13]]$median,
                           prediction_order[[j]][[14]]$median,
                           prediction_order[[j]][[15]]$median,
                           prediction_order[[j]][[16]]$median,
                           prediction_order[[j]][[17]]$median,
                           prediction_order[[j]][[18]]$median,
                           prediction_order[[j]][[19]]$median,
                           prediction_order[[j]][[20]]$median
                         ),
                         
                         x_coord = c(data_variables[[1]]$xcoord,
                                     data_variables[[2]]$xcoord,
                                     data_variables[[3]]$xcoord,
                                     data_variables[[4]]$xcoord,
                                     data_variables[[5]]$xcoord,
                                     data_variables[[6]]$xcoord,
                                     data_variables[[7]]$xcoord,
                                     data_variables[[8]]$xcoord,
                                     data_variables[[9]]$xcoord,
                                     data_variables[[10]]$xcoord,
                                     data_variables[[11]]$xcoord,
                                     data_variables[[12]]$xcoord,
                                     data_variables[[13]]$xcoord,
                                     data_variables[[14]]$xcoord,
                                     data_variables[[15]]$xcoord,
                                     data_variables[[16]]$xcoord,
                                     data_variables[[17]]$xcoord,
                                     data_variables[[18]]$xcoord,
                                     data_variables[[19]]$xcoord,
                                     data_variables[[20]]$xcoord),
                         
                         y_coord = c(data_variables[[1]]$ycoord,
                                     data_variables[[2]]$ycoord,
                                     data_variables[[3]]$ycoord,
                                     data_variables[[4]]$ycoord,
                                     data_variables[[5]]$ycoord,
                                     data_variables[[6]]$ycoord,
                                     data_variables[[7]]$ycoord,
                                     data_variables[[8]]$ycoord,
                                     data_variables[[9]]$ycoord,
                                     data_variables[[10]]$ycoord,
                                     data_variables[[11]]$ycoord,
                                     data_variables[[12]]$ycoord,
                                     data_variables[[13]]$ycoord,
                                     data_variables[[14]]$ycoord,
                                     data_variables[[15]]$ycoord,
                                     data_variables[[16]]$ycoord,
                                     data_variables[[17]]$ycoord,
                                     data_variables[[18]]$ycoord,
                                     data_variables[[19]]$ycoord,
                                     data_variables[[20]]$ycoord),
                         
                         bat = c(data_variables[[1]]$batimetria,
                                 data_variables[[2]]$batimetria,
                                 data_variables[[3]]$batimetria,
                                 data_variables[[4]]$batimetria,
                                 data_variables[[5]]$batimetria,
                                 data_variables[[6]]$batimetria,
                                 data_variables[[7]]$batimetria,
                                 data_variables[[8]]$batimetria,
                                 data_variables[[9]]$batimetria,
                                 data_variables[[10]]$batimetria,
                                 data_variables[[11]]$batimetria,
                                 data_variables[[12]]$batimetria,
                                 data_variables[[13]]$batimetria,
                                 data_variables[[14]]$batimetria,
                                 data_variables[[15]]$batimetria,
                                 data_variables[[16]]$batimetria,
                                 data_variables[[17]]$batimetria,
                                 data_variables[[18]]$batimetria,
                                 data_variables[[19]]$batimetria,
                                 data_variables[[20]]$batimetria),
                         time = rep(1:k, each = 10000)) 
}

for(i in 1:50) {
  print(i)
  data_all_2[[i]]$ocurrences <- ifelse(data_all_2[[i]]$probability >= 0.5, 1, 0)
  data_all_2[[i]]$fit_ocurrences <- ifelse(data_all_2[[i]]$pred_probability >= cutoff[i], 1, 0)
}

plot_pred <- list()
for(i in 1:50) { 
  plot_pred[[i]] <- ggplot(data = data_all_2[[i]]) + theme_bw() +
    geom_tile(aes(x=x_coord, y=y_coord, fill=pred_probability)) +
    labs(fill="probability", x = "", y = "") +
    scale_fill_viridis_c(option = "turbo",
                         limits = c(0,1)) + facet_wrap(~ time)
  }


confusion_all <- list()
# Calcular medidas
for (j in 1:50) {
  print(j)
  confusion_matrices <- list()
  for (year in unique(data_all_2[[j]]$time)) {
    # Filtrar los datos para el año actual
    subset_data <- data_all_2[[j]][data_all_2[[j]]$time == year, ]
    
    # Calcular la matriz de confusión para el año actual
    conf_matrix <- table(subset_data$ocurrences, subset_data$fit_ocurrences)
    
    # Almacenar la matriz de confusión en la lista
    confusion_matrices[[as.character(year)]] <- conf_matrix
  }
  # Almacenar la última matriz de confusión en la lista confusion_all
  confusion_all[[j]] <- confusion_matrices
}

sensitivity_all <- list()
specificity_all <- list()
accuracy_all <- list()

for(j in 1:50) {
  print(j)
  sensitivity <- list()
  specificity <- list()  
  accuracy <- list()
  
  tryCatch({
    for(i in 1:k) {
      sensitivity[[i]] <- sensitivity(confusion_all[[j]][[i]])
      specificity[[i]] <- specificity(confusion_all[[j]][[i]])
      accuracy[[i]] <- accuracy(data_all_2[[j]][data_all_2[[j]]$time == i,]$ocurrences,
                                data_all_2[[j]][data_all_2[[j]]$time == i,]$fit_ocurrences)
    }
    sensitivity_all[[j]] <- sensitivity
    specificity_all[[j]] <- specificity
    accuracy_all[[j]] <- accuracy
  }, error = function(e) {
    cat("Error in iteration", j, ":", conditionMessage(e), "\n")
    # You can choose to handle the error in a specific way or simply continue to the next iteration
  })
}


data_measures <- list()
for(i in 1:50) {
  print(i)
  tryCatch({
    data_measures[[i]] <- data.frame(
      sensitivity = unlist(sensitivity_all[[i]]),
      specificity = unlist(specificity_all[[i]]),
      accuracy = unlist(accuracy_all[[i]]),
      time = 1:20
    )
  }, error = function(e) {
    cat("Error in iteration", i, ":", conditionMessage(e), "\n")
    # You can choose to handle the error in a specific way or simply continue to the next iteration
  })
}

# Crea una lista con las medianas de cada columna
medians_list <- lapply(seq_len(20), function(i) {
  do.call(rbind, lapply(data_measures, function(df) df[i, ]))
})

# Combina la lista de medianas en un único data.frame
median_data_frame <- data.frame(
  sensitivity_median = sapply(medians_list, function(x) median(x$sensitivity)),
  specificity_median = sapply(medians_list, function(x) median(x$specificity)),
  accuracy_median = sapply(medians_list, function(x) median(x$accuracy)))

# Crea una lista con los cuantiles 2.5 de cada columna
q2.5_list <- lapply(seq_len(20), function(i) {
  do.call(rbind, lapply(data_measures, function(df) df[i, ]))
})

# Combina la lista de cuantiles 2.5 en un único data.frame
q2.5_data_frame <- data.frame(
  sensitivity_q2.5 = sapply(q2.5_list, function(x) quantile(x$sensitivity, prob = 0.025)),
  specificity_q2.5 = sapply(q2.5_list, function(x) quantile(x$specificity, prob = 0.025)),
  accuracy_q2.5 = sapply(q2.5_list, function(x) quantile(x$accuracy, prob = 0.025))
)

# Crea una lista con los cuantiles 97.5 de cada columna
q97.5_list <- lapply(seq_len(20), function(i) {
  do.call(rbind, lapply(data_measures, function(df) df[i, ]))
})

# Combina la lista de cuantiles 97.5 en un único data.frame
q97.5_data_frame <- data.frame(
  sensitivity_q97.5 = sapply(q97.5_list, function(x) quantile(x$sensitivity, prob = 0.975)),
  specificity_q97.5 = sapply(q97.5_list, function(x) quantile(x$specificity, prob = 0.975)),
  accuracy_q97.5 = sapply(q97.5_list, function(x) quantile(x$accuracy, prob = 0.975))
)

data_measures_all <- cbind(median_data_frame,
                           q2.5_data_frame,
                           q97.5_data_frame)
data_measures_all$time <- 1:20

# Transforma los datos al formato largo
data_measures_long <- tidyr::gather(data_measures_all,
                                    key = "Variable",
                                    value = "Value",
                                    sensitivity_median, specificity_median, accuracy_median,
                                    sensitivity_q2.5, specificity_q2.5, accuracy_q2.5,
                                    sensitivity_q97.5, specificity_q97.5, accuracy_q97.5)


# Suponiendo que ya tienes la columna "Tiempo" en tu data.frame
# Puedes ajustar esto según tus nombres de columna reales

# Crear el gráfico
ggplot(data_measures_long, aes(x = time, y = Value, color = Variable)) +
  geom_line(aes(group = interaction(Variable, ifelse(grepl("median", Variable), "median", ifelse(grepl("q2.5", Variable), "q2.5", "q97.5"))))) +
  geom_point(aes(shape = Variable), size = 2) +
  labs(x = "Tiempo", y = "Valor", color = "Variable", shape = "Variable") +
  scale_x_continuous(breaks = seq(0, 20, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_minimal()


# Crear el gráfico con "Time" en lugar de "Tiempo"
# Crear una columna para indicar el tipo de medida (accuracy, sensitivity, specificity)
data_measures_long$Measure_Type <- ifelse(grepl("accuracy", data_measures_long$Variable), "Accuracy",
                                          ifelse(grepl("sensitivity", data_measures_long$Variable), "Sensitivity", "Specificity"))

measures_plot_f <- ggplot(data_measures_long, aes(x = time, y = Value, color = Measure_Type, linetype = ifelse(grepl("median", Variable), "Median", ifelse(grepl("q2.5", Variable), "Q2.5", "Q97.5")))) +
  geom_line() +
  labs(x = "Time", y = "", color = "Measure", linetype = "Statistic") +
  scale_x_continuous(breaks = seq(0, 20, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("Accuracy" = "#1f78b4", "Sensitivity" = "#e41a1c", "Specificity" = "#33a02c")) +
  scale_linetype_manual(values = c("Median" = "solid", "Q2.5" = "dashed", "Q97.5" = "dashed")) +
  theme_minimal()


save.image(file = "./final_analysis.RData")

