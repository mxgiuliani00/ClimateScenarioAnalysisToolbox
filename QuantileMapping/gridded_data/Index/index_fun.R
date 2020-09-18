INDEX_data <- function(var, init_n, model, periods) {
  
  switch (var,
          'pav' = { A  <-  'pr';
          A1 <- 'Precipitation';
          u = 'mm/day'
          }, 'tas' = {A <- 'tas';
          A1 <- 'Temperature';
          u = 'K'
          },
          stop('define the variable')
  )
  
  wd <- paste(init_n, var, sep = '/')
  setwd(wd)
  
  y <- list()
  k <- 1
  dur <- 0
  for (p in periods) {
    setwd(paste(wd, p, sep = '/'))
    
    f <- list.files(pattern = '*.nc')
    name <- grepl(model, f)
    
    nc <- nc_open(f[name])
    x <- ncvar_get(nc, A)
    x <- apply(x,
               1,
               function(x) mean(x, na.rm = T)
    )
    
    nY <- length(x)/365
    dur <- dur + length(x)/365
    y[[k]] <- matrix(x, 365, nY)
    k <- k + 1
  }
  y <- do.call(cbind, y)
  
  yp <- list(x = y,
             Y = 2100 - dur + 1,
             datalab = paste(A1, u),
             titolo = paste(A1, 'MASH')
  )
  return(yp)
}

plot_index <- function(x, titolo, pss, model){
  
  nR <- nrow(x)
  nC <- ncol(x)
  
  
  if (nC > 4) {
    xm <- x[, c(1, ncol(x))]
  }else{
    xm <- x
  }
  
  xm <- cbind(
    apply(xm, 2, as.numeric), 
    rep(1, nrow(xm))
  )
  
  xm <- as.data.frame(xm)
  c3 <- paste(pss, model, sep = '_')
  colnames(xm) <- c('Years', 'index', 'V3')
  
  
  ggplot(data = xm, aes(Years, V3))+
    geom_raster(aes(fill = index))+
    scale_fill_viridis_c()+
    ylab(c3)+
    ggtitle(titolo)+
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
  ggsave(paste0("Index_",titolo, ".pdf"),
         plot = last_plot(), device = "pdf", scale = 1, dpi = 300)

}