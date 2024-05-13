##################################################
# BART protocol for GBIF and ISIMIP data
# Date: 29/02/2024
##################################################

##################################################
# 0. Load libraries, functions and supp. data
##################################################

# Load packages
library(tidyverse)
library(terra)
library(tidyterra)

# Load functions
source("R/invert_polygon.R")
source("R/download_data.R")
source("R/clean_coordinates.R")
source("R/generate_pseudo_absences.R")
source("R/prediction.R")
source("R/functional_responses.R")
source("R/cross_validation.R")

# Load world map
sf::sf_use_s2(FALSE)
global_land_mask <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") %>%
  sf::st_as_sfc() %>%
  sf::st_union() %>%
  sf::st_make_valid() %>%
  sf::st_wrap_dateline()

global_ocean_mask <- invert_polygon(global_land_mask)

##################################################
# 1. Download GBIF (presence-only) data
##################################################

# Config
user <- "albaf13" # your gbif.org username 
pwd <- "Minino2020" # your gbif.org password
email <-  "alba.fuster1398@gmail.com" # your email

sp_names <- c("Natator depressus",
              "Dermochelys coriacea",
              "Caretta caretta",
              "Lepidochelys olivacea",
              "Chelonia mydas",
              "Lepidochelys kempii",
              "Eretmochelys imbricata")

# Download GBIF data
occurrence_list <- list()
occurrence_list$raw_occur <- download_gbif_occurrences(sp_names,
                                                       user, pwd, email,
                                                       start_year = 2015,
                                                       end_year = 2024,
                                                       process_data = TRUE)

occurrence_list$raw_occur_pt <- download_gbif_occurrences(sp_names,
                                                       user, pwd, email,
                                                       start_year = 2015,
                                                       end_year = 2024,
                                                       process_data = TRUE)

##################################################
# 2. Clean coordinates
##################################################

# Remove NA coordinates
occurrence_list$noNA_occur <- lapply(occurrence_list$raw_occur, function(x) {
  x %>%
    dplyr::filter(!is.na(decimalLongitude)) %>%
    dplyr::filter(!is.na(decimalLatitude))
})

# Round coordinates to 4 digits
occurrence_list$round_occur <- lapply(occurrence_list$noNA_occur, function(x) {
  x %>% mutate(decimalLongitude = round(decimalLongitude, 4),
               decimalLatitude = round(decimalLatitude, 4))
})

# Remove duplicated points according to longitude and latitude (WGS84)
occurrence_list$uniq_occur <- lapply(occurrence_list$round_occur, function(x) {
  remove_duplicate_points(x, coords = c("decimalLongitude", "decimalLatitude"))
})

# Remove points outside the ocean (given a polygon of the land just keep points outside the polygon)
occurrence_list$sea_occur <- lapply(occurrence_list$uniq_occur, function(x) {
  remove_points_poly(x, global_land_mask, TRUE, c("decimalLongitude", "decimalLatitude")) 
})

# Filter with CoordinateCleaner tests
occurrence_list$clean_occur <- lapply(occurrence_list$sea_occur, function(x) {
  x <- CoordinateCleaner::clean_coordinates(x = x,
                         lon = "decimalLongitude",
                         lat = "decimalLatitude",
                         species ="species",
                         tests = c("institutions", "validity", "zeros"))
  x <- x[x$.summary, c("decimalLongitude", "decimalLatitude", "year", "species")]
  x
})

# Plot keeped and deleted points
for (sp in sp_names){
  plot(global_ocean_mask, col = "aliceblue", main = sp)
  plot(sf::st_geometry(sf::st_as_sf(occurrence_list$raw_occur[[sp]], coords = c("decimalLongitude", "decimalLatitude"))), col = "red", add = T)
  plot(sf::st_geometry(sf::st_as_sf(occurrence_list$clean_occur[[sp]], coords = c("decimalLongitude", "decimalLatitude"))), col = "green", add = T)
  
}

##################################################
# 3. Load environmental variables
##################################################

# List raster paths and file names
historical_dir <- "covariates/historical/"
projection_1.6_dir <- "covariates/projection_1.6/"
projection_5.8_dir <- "covariates/projection_5.8/"
historical_vars <- list.dirs(path = historical_dir)[-1]
projected_vars_1.6 <- list.dirs(path = projection_1.6_dir)[-1]
projected_vars_5.8 <- list.dirs(path = projection_5.8_dir)[-1]

env_var_names <- basename(historical_vars)

# Validate input rasters
stopifnot(basename(historical_vars) == basename(projected_vars_1.6)) # Check for same variables
stopifnot(basename(historical_vars) == basename(projected_vars_5.8)) # Check for same variables

years <- sapply(env_var_names, function(name) length(list.files(file.path(historical_dir, name))))
stopifnot(length(unique(years)) == 1) # Check for same number of files historical

years <- sapply(env_var_names, function(name) length(list.files(file.path(projection_1.6_dir, name))))
stopifnot(length(unique(years)) == 1) # Check for same number of files projections

years <- sapply(env_var_names, function(name) length(list.files(file.path(projection_5.8_dir, name))))
stopifnot(length(unique(years)) == 1) # Check for same number of files projections

# Load rasters from past (historical) and future (projections)
covariate_list <- list()
covariate_list$historical <- list()
for (i in seq_along(env_var_names)){
  covariate_list$historical[[env_var_names[i]]] <- terra::rast(list.files(file.path(historical_dir, env_var_names[i]), full.names = TRUE))
}

covariate_list$projected_1.6 <- list()
for (i in seq_along(env_var_names)){
  covariate_list$projected_1.6[[env_var_names[i]]] <- terra::rast(list.files(file.path(projection_1.6_dir, env_var_names[i]), full.names = TRUE))
}

covariate_list$projected_5.8 <- list()
for (i in seq_along(env_var_names)){
  covariate_list$projected_5.8[[env_var_names[i]]] <- terra::rast(list.files(file.path(projection_5.8_dir, env_var_names[i]), full.names = TRUE))
}

##################################################
# 4. Process covariates
##################################################

# Crop to globe extent and keep ocean values
for (i in seq_along(env_var_names)){
  # Crop to globe extent
  covariate_list$historical[[env_var_names[i]]] <- terra::crop(covariate_list$historical[[env_var_names[i]]], ext(-180, 180, -90, 90))
  covariate_list$projected_1.6[[env_var_names[i]]] <- terra::crop(covariate_list$projected_1.6[[env_var_names[i]]], ext(-180, 180, -90, 90))
  covariate_list$projected_5.8[[env_var_names[i]]] <- terra::crop(covariate_list$projected_5.8[[env_var_names[i]]], ext(-180, 180, -90, 90))
  # Extend to globe extent
  covariate_list$historical[[env_var_names[i]]] <- terra::extend(covariate_list$historical[[env_var_names[i]]], ext(-180, 180, -90, 90))
  covariate_list$projected_1.6[[env_var_names[i]]] <- terra::extend(covariate_list$projected_1.6[[env_var_names[i]]], ext(-180, 180, -90, 90))
  covariate_list$projected_5.8[[env_var_names[i]]] <- terra::extend(covariate_list$projected_5.8[[env_var_names[i]]], ext(-180, 180, -90, 90))
  # Remove land values
  covariate_list$historical[[env_var_names[i]]] <- terra::mask(covariate_list$historical[[env_var_names[i]]], terra::vect(global_ocean_mask))
  covariate_list$projected_1.6[[env_var_names[i]]] <- terra::mask(covariate_list$projected_1.6[[env_var_names[i]]], terra::vect(global_ocean_mask))
  covariate_list$projected_5.8[[env_var_names[i]]] <- terra::mask(covariate_list$projected_5.8[[env_var_names[i]]], terra::vect(global_ocean_mask))
}

# Scale rasters for model fitting
covariate_list$scl_historical <- list()
covariate_list$scl_projected_1.6 <- list()
covariate_list$scl_projected_5.8 <- list()
covariate_list$m_scl_historical <- list()

for (i in seq_along(env_var_names)){
  # Compute mean and sd of the historical time series for each environmental variable
  m_hist <- mean(as.vector(covariate_list$historical[[env_var_names[i]]]), na.rm = TRUE)
  sd_hist <- sd(as.vector(covariate_list$historical[[env_var_names[i]]]), na.rm = TRUE)
  # Scale historical and projected with historical mean and sd
  covariate_list$scl_historical[[env_var_names[i]]] <- terra::scale(covariate_list$historical[[env_var_names[i]]], center = m_hist, scale = sd_hist)
  covariate_list$scl_projected_1.6[[env_var_names[i]]] <- terra::scale(covariate_list$projected_1.6[[env_var_names[i]]], center = m_hist, scale = sd_hist)
  covariate_list$scl_projected_5.8[[env_var_names[i]]] <- terra::scale(covariate_list$projected_5.8[[env_var_names[i]]], center = m_hist, scale = sd_hist)
  # Get mean of scaled historical rasters
  covariate_list$m_scl_historical[[env_var_names[i]]] <- terra::mean(covariate_list$scl_historical[[env_var_names[i]]], na.rm = TRUE)
}

# Compute mean of historical variables (not scaled) for functional responses calculation
covariate_list$m_historical <- list()
for (i in seq_along(env_var_names)){
  covariate_list$m_historical[[env_var_names[i]]] <- terra::mean(covariate_list$historical[[env_var_names[i]]], na.rm = TRUE)
}

# Include longitude and latitude for native ranges modeling
data_latlong <- as.data.frame(covariate_list$m_scl_historical[[1]], xy=TRUE)[,-3]
data_latlong$decimalLongitude_scale <- data_latlong[,1]
data_latlong$decimalLatitude_scale <- data_latlong[,2]
raster_long <- terra::rast(data_latlong[,-4])
raster_long <- terra::extend(raster_long, covariate_list$m_scl_historical[[1]])
terra::crs(raster_long) <- "epsg:4326"
raster_long <- terra::scale(raster_long)
raster_lat <- terra::rast(data_latlong[,-3])
raster_lat <- terra::extend(raster_lat, covariate_list$m_scl_historical[[1]])
terra::crs(raster_lat) <- "epsg:4326"
raster_lat <- terra::scale(raster_lat)
covariate_list$m_scl_historical_nr <- c(terra::rast(covariate_list$m_scl_historical), raster_long, raster_lat) 

##################################################
# 5. Create model matrix and pseudoabsences
##################################################

# Create model matrix with occurrences and remove points with NA values in any environmental variable
model_data <- list()
model_data$fit_data <- lapply(occurrence_list$clean_occur, function(x) {
  fit_data <- as.data.frame(x[, c("decimalLongitude", "decimalLatitude", "year", "species")])
  design_matrix <- terra::extract(
    terra::rast(covariate_list$m_scl_historical),
    fit_data[, c("decimalLongitude", "decimalLatitude")]
  )
  
  fit_data <- cbind(fit_data, design_matrix) %>%
    tidyr::drop_na()
  fit_data
})

# Generate balanced random pseudoabsences
set.seed(4567)
model_data$fit_data_complete <- lapply(model_data$fit_data, function(x) {
  generate_pseudo_absences(x, global_ocean_mask, terra::rast(covariate_list$m_scl_historical))
})

# Non scale variables
model_data$fit_data_non_scale <- lapply(model_data$fit_data_complete, function(x){
  fit_data <- as.data.frame(x[, c("decimalLongitude", "decimalLatitude", "year", "species","pa")])
  design_matrix <- terra::extract(
    terra::rast(covariate_list$m_historical),
    fit_data[, c("decimalLongitude", "decimalLatitude")]
  )
  
  fit_data <- cbind(fit_data, design_matrix) %>%
    tidyr::drop_na()
  fit_data
})

# Data for native range
model_data$fit_data_nr <- lapply(model_data$fit_data_complete, function(x){
  fit_data <- as.data.frame(x[, c("decimalLongitude", "decimalLatitude", "year", "species","pa")])
  design_matrix <- terra::extract(
    covariate_list$m_scl_historical_nr,
    fit_data[, c("decimalLongitude", "decimalLatitude")]
  )
  
  fit_data <- cbind(fit_data, design_matrix) %>%
    tidyr::drop_na()
  fit_data
})

# Summary of model data
model_data_summary <- data.frame(species = as.character(), noNA = as.numeric(),
                                 unique = as.numeric(), sea = as.numeric(),
                                 presences = as.numeric(), pseudoabsences = as.numeric())
for (i in seq_along(sp_names)){
  sp <- sp_names[i]
  model_data_summary[i,] <- c(sp,
                              nrow(occurrence_list$noNA_occur[[sp]]), 
                              nrow(occurrence_list$uniq_occur[[sp]]),
                              nrow(occurrence_list$sea_occur[[sp]]),
                              sum(model_data$fit_data_complete[[sp]]$pa == 1),
                              sum(model_data$fit_data_complete[[sp]]$pa == 0))
  
  plot(global_ocean_mask, col = "aliceblue", main = sp)
  plot(sf::st_geometry(sf::st_as_sf(model_data$fit_data_complete[[sp]][model_data$fit_data_complete[[sp]]$pa == 0, ], coords = c("decimalLongitude", "decimalLatitude"))), col = "red", add = T)
  plot(sf::st_geometry(sf::st_as_sf(model_data$fit_data_complete[[sp]][model_data$fit_data_complete[[sp]]$pa == 1, ], coords = c("decimalLongitude", "decimalLatitude"))), col = "green", add = T)
  
}

##################################################
# 6. Fit the model
##################################################

# Spatial historical model ------------------------------------------------------- 
## Estimation
### Suitable habitat (scale)
set.seed(14523)
bart_model_sh_scl <- list()
for(i in 1:length(sp_names)) {
  print(i)
  bart_model_sh_scl[[i]] <- dbarts::bart(x.train = model_data$fit_data_complete[[i]][,env_var_names],
                                  y.train = model_data$fit_data_complete[[i]][,"pa"],
                                                keeptrees = TRUE)
}

for(i in 1:length(sp_names)) {
  print(i)
  invisible(bart_model_sh_scl[[i]]$fit$state) 
}

### Suitable habitat (non scale)
set.seed(14523)
bart_model_sh <- list()
for(i in 1:length(sp_names)) {
  print(i)
  bart_model_sh[[i]] <- dbarts::bart(x.train = model_data$fit_data_non_scale[[i]][,env_var_names],
                                       y.train = model_data$fit_data_non_scale[[i]][,"pa"],
                                       keeptrees = TRUE)
}

for(i in 1:length(sp_names)) {
  print(i)
  invisible(bart_model_sh[[i]]$fit$state) 
}

### Native range (scale)
set.seed(14523)
bart_model_nr <- list()
for(i in 1:length(sp_names)) {
  print(i)
  bart_model_nr[[i]] <- dbarts::bart(x.train = model_data$fit_data_nr[[i]][,c(env_var_names,
                                                                                    "decimalLongitude_scale",
                                                                                    "decimalLatitude_scale")],
                                     y.train = model_data$fit_data_nr[[i]][,"pa"],
                                     keeptrees = TRUE)
}

for(i in 1:length(sp_names)) {
  print(i)
  invisible(bart_model_nr[[i]]$fit$state) 
}

## Prediction
### Suitable habitat (scale)
bart_pred_sh_scl <- list()
for(i in 1:length(sp_names)) {
  print(i)
  bart_pred_sh_scl[[i]] <- prediction_BART(bart_model = bart_model_sh_scl[[i]],
                                    stack_cov = terra::rast(covariate_list$m_scl_historical))
}

for(i in 1:length(sp_names)) {
  print(i)
  terra::writeRaster(bart_pred_sh_scl[[i]],
                     paste0("./results/IPSL/historical/suitable_habitat/",sp_names[i],".tiff")) 
}

p_pred_sh_scl <- list()
for(i in 1:length(sp_names)) {
  p_pred_sh_scl[[i]] <- ggplot() +
  geom_spatraster(data = bart_pred_sh_scl[[i]]) +
  facet_wrap(~lyr) + scale_fill_gradientn(colours = c("#A1D4B1","#2BAF90",
                                                      "#F1A512","#DD4111",
                                                      "#8C0027"),
                                       limits = c(0,1), name = "Probability") +
  geom_sf(data = global_land_mask, fill = "antiquewhite") +
  theme_minimal()
}

### Native ranges (scale)
bart_pred_nr_scl <- list()
for(i in 1:length(sp_names)) {
  print(i)
  bart_pred_nr_scl[[i]] <- prediction_BART(bart_model = bart_model_nr[[i]],
                                           stack_cov = covariate_list$m_scl_historical_nr)
}

for(i in 1:length(sp_names)) {
  print(i)
  terra::writeRaster(bart_pred_nr_scl[[i]],
                     paste0("./results/IPSL/historical/native_range/",sp_names[i],".tiff")) 
}

p_pred_nv_scl <- list()
for(i in 1:length(sp_names)) {
  p_pred_nv_scl[[i]] <- ggplot() +
    geom_spatraster(data = bart_pred_nr_scl[[i]]) +
    facet_wrap(~lyr) + scale_fill_gradientn(colours = c("#A1D4B1","#2BAF90",
                                                        "#F1A512","#DD4111",
                                                        "#8C0027"),
                                            limits = c(0,1), name = "Probability") +
    geom_sf(data = global_land_mask, fill = "antiquewhite") +
    theme_minimal()
}

## Projections (historical and future)  -------------------------------------------------------
### Stack with the covariates order by year from 1950 to 2014
tmp_stack_historical <- list()
for(i in 1:dim(covariate_list$scl_historical[[1]])[3]) {
  tmp_stack_historical[[i]] <- terra::rast(lapply(covariate_list$scl_historical, function(x) {
    x[[i]]
  }))
}

### prediction historical
bart_proj_hist_all <- list()
bart_proj_hist <- list()
for(i in 1:length(sp_names)) {
  print(i)
  for(j in 1:dim(covariate_list$scl_historical[[1]])[3]) {
    print(j)
    bart_proj_hist[[j]] <- prediction_BART(bart_model = bart_model_sh_scl[[i]],
                                      stack_cov = tmp_stack_historical[[j]])
  }
  bart_proj_hist_all[[i]] <- bart_proj_hist
}

bart_proj_hist_all_species <- list()
for(i in 1:length(sp_names)) {
  print(i)
  bart_proj_hist_all_species[[i]] <- terra::rast(bart_proj_hist_all[[i]])
}

for(i in 1:length(sp_names)) {
  print(i)
  terra::writeRaster(bart_proj_hist_all_species[[i]],
                     paste0("./results/IPSL/historical/projections/",sp_names[i],".tiff")) 
}

p_pred_hist_scl <- list()
p_pred_hist_all_scl <- list()
for(i in 1:length(sp_names)) {
  print(i)
  for(j in 1:dim(covariate_list$scl_historical[[1]])[3]) {
    print(j)
    p_pred_hist_scl[[j]] <- ggplot() +
      geom_spatraster(data = bart_proj_hist_all[[i]][[j]]) +
      facet_wrap(~lyr) + scale_fill_gradientn(colours = c("#A1D4B1","#2BAF90",
                                                          "#F1A512","#DD4111",
                                                          "#8C0027"),
                                              limits = c(0,1), name = "Probability") +
      geom_sf(data = world, fill = "antiquewhite") +
      theme_minimal() 
  }
  p_pred_hist_all_scl[[i]] <- p_pred_hist_scl
}

### Stack with the covariates order by year from 2015 to 2100
tmp_stack_future_1.6 <- list()
for(i in 1:dim(covariate_list$scl_projected_1.6[[1]])[3]) {
  tmp_stack_future_1.6[[i]] <- terra::rast(lapply(covariate_list$scl_projected_1.6, function(x) {
    x[[i]]
  }))
}

tmp_stack_future_5.8 <- list()
for(i in 1:dim(covariate_list$scl_projected_5.8[[1]])[3]) {
  tmp_stack_future_5.8[[i]] <- terra::rast(lapply(covariate_list$scl_projected_5.8, function(x) {
    x[[i]]
  }))
}


### prediction future 1.6
bart_proj_all_1.6 <- list()
bart_proj_1.6 <- list()
for(i in 1:length(sp_names)) {
  print(i)
  for(j in 1:dim(covariate_list$scl_projected_1.6[[1]])[3]) {
    print(j)
    bart_proj_1.6[[j]] <- prediction_BART(bart_model = bart_model_sh_scl[[i]],
                                      stack_cov = tmp_stack_future_1.6[[j]])
  }
  bart_proj_all_1.6[[i]] <- bart_proj_1.6
}

bart_proj_1.6_all_species <- list()
for(i in 1:length(sp_names)) {
  print(i)
  bart_proj_1.6_all_species[[i]] <- terra::rast(bart_proj_all_1.6[[i]])
}

for(i in 1:length(sp_names)) {
  print(i)
  terra::writeRaster(bart_proj_1.6_all_species[[i]],
                     paste0("./results/IPSL/projected_1.6/",sp_names[i],".tiff")) 
}

p_pred_1.6_scl <- list()
p_pred_1.6_all_scl <- list()
for(i in 1:length(sp_names)) {
  print(i)
  for(j in 1:dim(covariate_list$scl_projected_1.6[[1]])[3]) {
    print(j)
    p_pred_1.6_scl[[j]] <- ggplot() +
      geom_spatraster(data = bart_proj_all_1.6[[i]][[j]]) +
      facet_wrap(~lyr) + scale_fill_gradientn(colours = c("#A1D4B1","#2BAF90",
                                                          "#F1A512","#DD4111",
                                                          "#8C0027"),
                                              limits = c(0,1), name = "Probability") +
      geom_sf(data = world, fill = "antiquewhite") +
      theme_minimal() 
  }
  p_pred_1.6_all_scl[[i]] <- p_pred_1.6_scl
}

### prediction future 5.8
bart_proj_all_5.8 <- list()
bart_proj_5.8 <- list()
for(i in 1:length(sp_names)) {
  print(i)
  for(j in 1:dim(covariate_list$scl_projected_5.8[[1]])[3]) {
    print(j)
    bart_proj_5.8[[j]] <- prediction_BART(bart_model = bart_model_sh_scl[[i]],
                                          stack_cov = tmp_stack_future_5.8[[j]])
  }
  bart_proj_all_5.8[[i]] <- bart_proj_5.8
}

bart_proj_5.8_all_species <- list()
for(i in 1:length(sp_names)) {
  print(i)
  bart_proj_5.8_all_species[[i]] <- terra::rast(bart_proj_all_5.8[[i]])
}

for(i in 1:length(sp_names)) {
  print(i)
  terra::writeRaster(bart_proj_5.8_all_species[[i]],
                     paste0("./results/IPSL/projected_5.8/",sp_names[i],".tiff")) 
}

p_pred_5.8_scl <- list()
p_pred_5.8_all_scl <- list()
for(i in 1:length(sp_names)) {
  print(i)
  for(j in 1:dim(covariate_list$scl_projected_5.8[[1]])[3]) {
    print(j)
    p_pred_5.8_scl[[j]] <- ggplot() +
      geom_spatraster(data = bart_proj_all_5.8[[i]][[j]]) +
      facet_wrap(~lyr) + scale_fill_gradientn(colours = c("#A1D4B1","#2BAF90",
                                                          "#F1A512","#DD4111",
                                                          "#8C0027"),
                                              limits = c(0,1), name = "Probability") +
      geom_sf(data = world, fill = "antiquewhite") +
      theme_minimal() 
  }
  p_pred_5.8_all_scl[[i]] <- p_pred_5.8_scl
}

# Functional responses -------------------------------------------------------
fr_turtles <- list()
for(i in 1:length(sp_names)) {
  print(i)
  fr_turtles[[i]] <- cov_relations(bart_model = bart_model_sh[[i]],
                              names = env_var_names,
                              data = model_data$fit_data_non_scale[[i]])
}

p <- list()
p_fr <- list()
for(j in 1:length(sp_names)) {
  print(j)
  for (i in 1:length(env_var_names)) {
    print(i)
   p[[i]] <- ggplot(fr_turtles[[j]][[i]], aes(x = var, y = prob)) +
      geom_line() +
      geom_line(aes(y = q25), linetype = "dashed") +
      geom_line(aes(y = q975), linetype = "dashed") +
      xlab(env_var_names[i]) +
      ylab("Probability") + theme_minimal()
  }
  p_fr[[j]] <- p 
}

# Cross validation -------------------------------------------------------
cv_turtles <- list()
for(i in 1:length(sp_names)) {
  print(i)
  cv_turtles[[i]] <- cv_bart(pa_coords = model_data$fit_data_complete[[i]],
                             layers = terra::rast(covariate_list$m_scl_historical), k = 10, seed = 16332)
}


# Calculate cutoff -------------------------------------------------------
pred <- list()
pred_mean <- list()
pred_loc <- list()
for(i in 1:length(sp_names)) {
  print(i)
  pred[[i]] <- dbarts:::predict.bart(bart_model_sh_scl[[i]],
                                     newdata = model_data$fit_data_complete[[i]][,env_var_names])
  pred_mean[[i]] <- colMeans(pred[[i]])
  pred_loc[[i]] <- data.frame(long = model_data$fit_data_complete[[i]]$decimalLongitude,
                              lat = model_data$fit_data_complete[[i]]$decimalLatitude,
                              pa = model_data$fit_data_complete[[i]]$pa,
                              prob = pred_mean[[i]])
  
}

cutoff <- list()
for(i in 1:length(sp_names)) {
  cutoff[[i]] <- optimalCutoff(actuals = pred_loc[[i]]$pa,
                               predictedScores = pred_loc[[i]]$prob,
                               optimiseFor = "misclasserror", returnDiagnostics = FALSE)
}

## cutoff GFDL y IPSL
##[[1]]
##[1] 0.5381179
##
##[[2]]
##[1] 0.4774161
##
##[[3]]
##[1] 0.5685869
##
##[[4]]
##[1] 0.5267385
##
##[[5]]
##[1] 0.4878081
##
##[[6]]
##[1] 0.6886046
##
##[[7]]
##[1] 0.5988242


## Functional group 
historical_sh_GFDL <- list()
for(i in 1:length(sp_names)) {
  historical_sh_GFDL[[i]] <- terra::rast(
    paste0("./results/GFDL/historical/suitable_habitat/",sp_names[i],".tiff")
  )
}

historical_sh_GFDL_mean <- list()
for(i in 1:length(sp_names)) {
  historical_sh_GFDL_mean[[i]] <- historical_sh_GFDL[[i]]["mean"]
}

hist_sh_GFDL_mean_all <- terra::mean(terra::rast(historical_sh_GFDL_mean))

pred_loc_all <- do.call(rbind,pred_loc)
pred_loc_all <- pred_loc_all[,-4]
pred_loc_all_sort <- pred_loc_all[order(pred_loc_all$pa, decreasing = TRUE),]
pred_loc_all_sort_duplicate <- pred_loc_all_sort[!duplicated(pred_loc_all_sort[, c("long","lat")]), ]
probability_fg <- terra::extract(hist_sh_GFDL_mean_all$mean, cbind(pred_loc_all_sort_duplicate$long,pred_loc_all_sort_duplicate$lat))
pred_loc_all_sort_duplicate$prob <- probability_fg
fg_cutoff <- optimalCutoff(actuals = pred_loc_all_sort_duplicate$pa, predictedScores = pred_loc_all_sort_duplicate$prob)

historical_sh_IPSL <- list()
for(i in 1:length(sp_names)) {
  historical_sh_IPSL[[i]] <- terra::rast(
    paste0("./results/IPSL/historical/suitable_habitat/",sp_names[i],".tiff")
  )
}

historical_sh_IPSL_mean <- list()
for(i in 1:length(sp_names)) {
  historical_sh_IPSL_mean[[i]] <- historical_sh_IPSL[[i]]["mean"]
}

historical_sh_IPSL_mean_all <- terra::mean(terra::rast(historical_sh_IPSL_mean))

pred_loc_all <- do.call(rbind,pred_loc)
pred_loc_all <- pred_loc_all[,-4]
pred_loc_all_sort <- pred_loc_all[order(pred_loc_all$pa, decreasing = TRUE),]
pred_loc_all_sort_duplicate <- pred_loc_all_sort[!duplicated(pred_loc_all_sort[, c("long","lat")]), ]
probability_fg <- terra::extract(historical_sh_IPSL_mean_all$mean, cbind(pred_loc_all_sort_duplicate$long,pred_loc_all_sort_duplicate$lat))
pred_loc_all_sort_duplicate$prob <- probability_fg
fg_cutoff <- optimalCutoff(actuals = pred_loc_all_sort_duplicate$pa, predictedScores = pred_loc_all_sort_duplicate$prob)
# 0.1318681 # GFDL
# 0.1305587 # IPSL