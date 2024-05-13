#*********************************
# SPATIO-TEMPORAL SIMULATION     #
# Biomass and sampling           #
#*********************************
# Alba Fuster-Alonso             #
#*********************************

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 

# Packages ---------------------------
library(INLA)
library(fields)
library(lattice)
library(reshape)
library(gridExtra)
library(RColorBrewer)
library(raster)
library(ggplot2)

# Function INLA for spatio-temporal term ---------------------------
source("./function.R")

# Parameters ---------------------------
coord1 <- 0 # coordinates
coord2 <- 10 # coordinates
coord3 <- 0 # coordinates
coord4 <- 10 # coordinates
variance <- 1 # spatial effect variance
kappa <- 0.5 # range = sqrt(8)/kappa
rho_1 <- 0.1 # Rho (temporal correlation)
rho_2 <- 0.7 # Rho (autoregresivo)
k <- 20 # k indica el numero de anyos
beta_0 <- 1 # intercept
beta_1 <- 7.5 # coefficient (degree 1) for the covariate (bathymetry)
#beta_2 <- 2.1 # coefficient (degree 2) for the covariate (bathymetry)
beta_3 <- -0.8
set.seed(852345)
vector_tiempo <- as.vector(arima.sim(list(order=c(1,0,0), ar=rho_2), n=k))
m <- 50 # Number of sampling points

# Simulated study area ---------------------------
campo_estudio <- function(coord1, coord2, coord3, coord4) {
  xy <- expand.grid(
    seq(coord1, coord2, length.out = (100)),
    seq(coord3, coord4, length.out = (100))
  )
  x <- as.numeric(xy[, 1])
  y <- as.numeric(xy[, 2])
  loc_xy <- cbind(x, y)

  return(loc_xy)
}

set.seed(852345) 
loc_xy <- campo_estudio(coord1, coord2, coord3, coord4)

# Simulation of the predictor terms ---------------------------
simulacion_variables <- function(variance, kappa, rho, k) {
  # Mesh
  prmesh1 <- inla.mesh.2d(loc = loc_xy, max.edge = c(0.4, 1))
  
  # Parameters for spatial effect
  params <- c(variance = variance, kappa = kappa)
  
  # INLA function
  set.seed(852345) 
  x_k <- book.rspde(loc_xy,
                    range = sqrt(8) / params[2],
                    sigma = sqrt(params[1]), n = k, mesh = prmesh1,
                    return.attributes = TRUE,
                    seed = 2468 
  )
  # Add temporal correlation
  x <- x_k
  for (j in 2:k) {
    x[, j] <- rho * x[, j - 1] + sqrt(1 - rho^2) * x_k[, j]
  }

  # Function for bathymetry
  batimetria <- function(x, y) {
    x_bat1 <- log(x * y + 1) * 100
    return(x_bat1)
  }
  
  # Function for bathymetry
  x_bat1 <- batimetria(loc_xy[, 1], loc_xy[, 2])
  
  # Scale bathymetry
  x_bat <- scale(x_bat1)
  
  # Function for temperature
  temperature <- function(y) {
    x_temp1 <- sqrt(y + 1) + 10
    return(x_temp1)
  }

  x_temp1 <- list()
  for (i in 1:k) {
    if (i == 1) {
      x_temp1[[i]] <- as.vector(scale(temperature(loc_xy[,2]) + 0.5))
    } else {
      x_temp1[[i]] <- x_temp1[[i-1]] + 0.5
    }
  }
  
  data_temp <- list()
  for(i in 1:k) {
    data_temp[[i]] <- data.frame(x = loc_xy[,1], y = loc_xy[,2], temp = x_temp1[[i]])
  }
  
  # Inicializar la matriz
  x_temp <- matrix(NA, nrow = 10000, ncol = 20)
  
  # Llenar la matriz con los valores de x_temp1[[i]]
  for (i in 1:20) {
    # Verificar si x_temp1[[i]] es una lista
    if (is.list(x_temp1[[i]]) && "temp" %in% names(x_temp1[[i]])) {
      x_temp[, i] <- x_temp1[[i]]$temp
    } else {
      # Si x_temp1[[i]] no es una lista o no tiene un componente 'temp', ajusta el código según tu estructura
      x_temp[, i] <- x_temp1[[i]]
    }
  }
  
  # Data.frame with predictor's terms
  variables <- data.frame(x = x, x_bat1 = x_bat1, x_bat = x_bat,
                          x_temp = x_temp)
  return(variables)
}

set.seed(852345) 
variables <- simulacion_variables(variance, kappa, rho = rho_1, k)

# Variable response simulation ---------------------------
simulacion_respuesta <- function(beta_0, beta_1, beta_2, beta_3, n) {

  lin_pred_mu <- list()
  linear_combination <- list()
  for (i in 1:k) {
    print(i)
    linear_combination[[i]] <- beta_0 +
      beta_1 * variables$x_bat +
      beta_3 * variables[,i+22] +
      variables[,i] +
      vector_tiempo[i]
    
    lin_pred_mu[[i]] <- boot::inv.logit(linear_combination[[i]])
  }

  # Simulation of the biomass with a bernouilli distribution
  presences_real <- list()
  for (i in 1:k) {
    presences_real[[i]] <- rbinom(n, 1, lin_pred_mu[[i]])
  }

  # Datas with biomass, bathymetry and locations
  data_variables <- list()
  for (i in 1:k) {
    data_variables[[i]] <- data.frame(
      probability = lin_pred_mu[[i]],
      batimetria = variables$x_bat,
      temp = variables[,i+22],
      xcoord = loc_xy[, 1],
      ycoord = loc_xy[, 2],
      ocurrences = presences_real[[i]]
    )
  }
  
  return(data_variables)
}

n <- as.numeric(dim(variables)[1]) # number of rows of variables's data frame

set.seed(852345) 
data_variables <- simulacion_respuesta(beta_0, beta_1, beta_2,beta_3, n)

# RData --------------------------- 
save.image(file = "./simulated_data.RData")

# Sampling ---------------------------
## Random sampling ---------------------------
muestreo_random <- function(m, q_random) {
  
  # Number of points selected
  data_random <- list()
  for (i in 1:k) {
    # Filtrar índices de presencias (ocurrences == 0) y ausencias (ocurrences == 1)
    indices_presencias <- which(data_variables[[i]]$ocurrences == 0)
    indices_ausencias <- which(data_variables[[i]]$ocurrences == 1)
    
    # Seleccionar m/2 puntos aleatorios de presencias y m/2 de ausencias
    isel_presencias <- sample(indices_presencias, m/2)
    isel_ausencias <- sample(indices_ausencias, m/2)
    
    # Combinar los índices seleccionados
    isel <- c(isel_presencias, isel_ausencias)
    
    # Almacenar los datos seleccionados
    data_random[[i]] <- data_variables[[i]][isel, ]
  }
  
  # Final data frame
  data_random_final <- data.frame()
  
  for (i in 1:k) {
    data_random_final <- rbind(data_random_final, data_random[[i]])
  }
  
  # Add time 
  data_random_final$time <- rep(1:k, each = m)
  
  return(data_random_final)
}


data_random_final <- list()
qseed <- list()
for(i in 1:50) {
  print(i)
  qseed[[i]] <- sample(1:10000,1)
  set.seed(qseed[[i]])
  data_random_final[[i]] <- muestreo_random(m, q_random)
}

# RData --------------------------- 
save.image(file = "./simulated_data_sampling.RData")
