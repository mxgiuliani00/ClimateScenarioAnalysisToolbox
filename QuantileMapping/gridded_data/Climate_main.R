setwd("~/Desktop/Gridded")

library(ncdf4)
library(raster)
library(rgdal)
library(lubridate)
library(geosphere)
library(reshape2)
library(ggplot2)
library(dplyr)
source('Climate_fun.R')

# Downscaling----

# a- couple control period
init <- '/Users/alessandroamaranto/Desktop/Gridded/Data/Down_in'
vars <- c('pav', 'tas')
id <- c('rca4') # add here if more models are investigated
ps <- 'cp'

X_c <- matrix(list())
i <- 1
for (var in vars) {
  X_c[[i]] <- return_control(init, var, ps, id)
  i <- i + 1
}

# b- couple future period
pss <- c('rcp45') # add here if more scenarios are investigated
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
init_n <- '/Users/alessandroamaranto/Desktop/Gridded/Data/Down_out'
i <- 1
for (var in vars) {
  nc_from_ts(X[[i]], init, var, id, pss, init_n)
  i <- i + 1
}

# Postprocessing and visualization----

# a- many possible futures
labs <- c('control', pss)
i <- 1
p <- list()
for (var in vars) {
  p[[i]] <- futures_plot_data(init_n, var, labs, id)
  i <- i+1
}

plot_futures(p, vars, labs, id)

# b- MESHdfs
model <- 'rca4'
periods <- c('control', 'rcp45')

X <- list()
i <- 1
for (var in vars) {
  X[[i]] <- MASH_data(var, init_n, model, periods)
  MASH_plot(X[[i]]$x, X[[i]]$Y, 10, 10, X[[i]]$datalab, X[[i]]$titolo)
  i <- i + 1
}






