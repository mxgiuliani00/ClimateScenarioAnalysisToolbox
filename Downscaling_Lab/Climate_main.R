library(ncdf4)
library(raster)
library(rgdal)
library(lubridate)
library(geosphere)
library(reshape2)
library(ggplot2)
source('ecqc.R')
source('ecui.R')
source('rclimdex.R')
source('rclimdexpatch.R')
source('Climate_fun.R')

# Downscaling----

# a- couple control period
init <- '/Users/alessandroamaranto/Desktop/Downscaling_Lab/Data/Down_in'
vars <- c('pav', 'tas')
id <- c('racmo', 'rca4')
ps <- 'cp'

X_c <- matrix(list())
i <- 1
for (var in vars) {
  X_c[[i]] <- return_control(init, var, ps, id)
  i <- i + 1
}

# b- couple future period
pss <- c('rcp45', 'rcp85')
X_f <- list()
i <- 1
for (var in vars) {
  X_f[[i]] <- return_future(init, var, pss, id)
  i <- i + 1
}

# c- downscale
X <- list()
i <- 1
for (var in vars) {
  X[[i]] <- downscale(var, X_f[[i]], X_c[[i]], init)
  i <- i + 1
}

# d- convert to netcdf
i <- 1
for (var in vars) {
  nc_from_ts(X[[i]], init, var, id, pss)
  i <- i + 1
}

# Postprocessing and visualization----

# a- many possible futures
labs <- c('control', pss)
init_n <- '/Users/alessandroamaranto/Desktop/Downscaling_Lab/Data/Down_out'
i <- 1
p <- list()
for (var in vars) {
  p[[i]] <- futures_plot_data(init_n, var, labs, id)
  i <- i+1
}

plot_futures(p, vars, labs, id)

# b- MESHdfs
model <- 'rca4'
periods <- c('control', 'rcp85')

X <- list()
i <- 1
for (var in vars) {
  X[[i]] <- MASH_data(var, init_n, model, periods)
  MASH_plot(X[[i]]$x, X[[i]]$Y, 10, 10, X[[i]]$datalab, X[[i]]$titolo)
  i <- i + 1
}


# c- climate index
rclimdex.env <- new.env()
rclimdex.env$parameter$bp.first <- 1988
rclimdex.env$parameter$bp.last <- 2005
rclimdex.env$parameter$lat <- 0
rclimdex.env$parameter$station <- 'not specified'
rclimdex.env$parameter$file.loc <- 'tba'
rclimdex.env$parameter$cal.type <- 'gregorian'

j <-1
pss <- c('rcp45', 'rcp85')
t <- seq.Date(as.Date('1988-01-01'), as.Date('2100-12-31'), 'days')
d <- day(t)
m <- month(t)
idx <- !(d == 29 & m == 2)
t <- t[idx]

nD <- length(t) 
I <- list()

for (p in pss) {
  periods <- c('control', pss[j])
  Y <- matrix(NA, nD, 6)
  i <- 1;
  for (var in vars) {
    X <- MASH_data(var, init_n, model, periods)
    x <- X$x
    Y[, i+3] <- as.vector(x)
    i <- i + 1
  }
  Y[ ,6] <- Y[ ,5]
  Y[ ,1] <- year(t)
  Y[ ,2] <- month(t)
  Y[ ,3] <- day(t)
  
  colnames(Y) <- c('year', 'month', 'day', 'prcp', 'tmin', 'tmax')
  rclimdex.env$data <- Y
  Input <- ClimdexInput()
  I[[j]] <- IndicesCalculation(Input)
  j <- j + 1
}  





