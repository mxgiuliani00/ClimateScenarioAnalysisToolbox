setwd("~/Desktop/Gridded/Index")
library(lubridate)
library(ncdf4)
source('index_fun.R')
source('ecqc.R')
source('ecui.R')
source('rclimdex.R')
source('rclimdexpatch.R')
library(ggplot2)

# a- climate index
init_n <- '/Users/alessandroamaranto/Desktop/Gridded/Data/Down_out'
vars <- c('pav', 'tas')
pss <- c('rcp45') # add here for other scenarios
model <- 'rca4'

rclimdex.env <- new.env()
rclimdex.env$parameter$bp.first <- 1988
rclimdex.env$parameter$bp.last <- 2005
rclimdex.env$parameter$lat <- 0
rclimdex.env$parameter$station <- 'not specified'
rclimdex.env$parameter$file.loc <- 'tba'
rclimdex.env$parameter$cal.type <- 'gregorian'

t <- seq.Date(as.Date('1988-01-01'), as.Date('2100-12-31'), 'days')
d <- day(t)
m <- month(t)
idx <- !(d == 29 & m == 2)
t <- t[idx]

nD <- length(t) 
I <- list()

j <-1
for (p in pss) {
  periods <- c('control', pss[j])
  Y <- matrix(NA, nD, 6)
  i <- 1;
  for (var in vars) {
    X <- INDEX_data(var, init_n, model, periods)
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

# Plot the indices
titles <- c('tx90p', "wsdi", "csdi", "r95p")
k <- 1

for (title in titles) {
  
  x <- lapply(I, function(l) l[[k]])
  x <- x[[1]] # temporary

  plot_index(x, title, pss, model)
  k <- k+1
  
}





