# Downscale----

return_control <- function(init, var, ps, id) {
  
  # return_control <- function(init, var, ps, id)
  # extract from a CORDEX netcdf file the nonnulls cells (i.e. those within the basin), and couples them with the
  # closest point in an observation netcdf file
  # args: 
  # init = initial path to the folder where the data are stored
  # var = variable name [one among pav or tas]
  # ps = scenario abbreviation [one among cp, rcp26, rcp45 and rcp85]
  # id = models ids
  # returns:
  # X_C: a list (length = number of variables) of matrices. Each list includes as mani matrices
  # as the climate models. Each matrix is a TxS, where T is the number of time step, and S are the 
  # number of cells within the domain
  
  
  # Set working directory as the investiated variable
  fin <- paste(init, var, sep = '/')
  setwd(fin)
  
  # Open observation file and extract lat, lon, vars
  fname <- list.files(pattern = '*.nc')
  nc <- nc_open(fname)
  
  lon <- ncvar_get(nc, 'lon')
  lat <- ncvar_get(nc, 'lat')
  
  # Check is a cell is nonnull
  x <- ncvar_get(nc, var)
  idx <- apply(
    x,
    2,
    function(x) sum(
      !is.na(x)
    )
  )
  idx <- idx != 0
  
  # Extract only internal points
  st <- data.frame(lon[idx], lat[idx])
  
  # List cordex cp files (one per model)
  setwd(paste(fin, ps, sep = '/'))
  fname <- list.files(pattern = '*.nc')
  
  k <- 1
  idx <- matrix(NA, length(id), length(ps))
  
  # Check is data for a certain model are not available
  for (f in id) {
    name <- grepl(f, fname)
    
    check <- sum(name) > 0
    if (check) {
      idx[k] <- k
    }
    k <- k + 1
  }
  
  idx <- idx[!is.na(idx)]
  
  Xc <- matrix(list(), length(id), length(ps))
  
  switch (var,
          'pav' = { A  = 'pr'
          }, 'tas' = {A = 'tas'
          },
          stop('define the case')
  )
  
  # Iterate through models and associate to each CORDEX cell the closest observation cell
  q <- 1
  m <- 1
  for (f in fname) {
    nc <- nc_open(f)
    if (q == 1) {
      
      lc <- ncvar_get(nc, 'lon')
      lac <- ncvar_get(nc, 'lat')
      
      p <- data.frame(lc, lac)
      
      d_min <- apply(st,
                     1,
                     function(x) which.min(
                       distCosine(x, p)
                     )
      )
      q <- q+1
    }
    
    y <- ncvar_get(nc, A)
    x <- matrix(NA, nrow(y), nrow(st))
    
    for (j in 1:ncol(x)) {
      x[ ,j] <- y[, d_min[j]]
    }
    
    Xc[[idx[m]]] <- x
    m <- m + 1
    nc_close(nc)
  }
  return(Xc)
}

return_future <- function(init, var, pss, id) {
  
  # return_future <- function(init, var, pss, id) 
  # function to couple each cell in each scenario (and model) with the closest cell in a netcdf of observations
  # args:
  #
  # init = initial path to the folder where the data are stored
  # var = variable name [one among pav or tas]
  # ps = scenario abbreviation [one among cp, rcp26, rcp45 and rcp85]
  # id = models ids
  # Returns:
  #
  # X_C: a list (length = number of variables) of matrices. Each list includes as mani matrices
  # as the climate models x the number of scenarios. Each matrix is a TxS, where T is the number of time step, and S are the 
  # number of cells within the domain
  
  Xv <- matrix(list(), length(pss), length(id))
  m <- 1
  
  for (ps in pss) {
    x <- return_control(init, var, ps, id)
    j <- 1
    while(j <= length(x)) {
      
      if (!is.null(x[[j]])) {
        Xv[[m, j]] <- x[[j]]
      }
      j <- j + 1
    }
    m <- m + 1
  }
  return(Xv)
}

downscale <- function(var, Xf, Xc, init) {
  # downscale <- function(var, Xf, Xc, init) 
  
  # Input:
  # variable name [pav or tas]
  # Xf = list of matrixes containing cordex projections
  # Xc = list of matrixes containing control period data
  # Return: 
  # A list of matrixes containing downscaled data
  
  switch (var,
          'pav' = { A <- 'pav';
          A1 <- 'pr';
          A2 <- 'P';
          mult <- 86400
          }, 'tas' = {A <- 'tas';
          A1 <- 'tas';
          A2 <- 'T';
          mult <- 1
          },
          stop('define the case')
  )
  
  fin <- paste(init, var, sep = '/')
  setwd(fin)
  
  fname <- list.files(pattern = '*.nc')
  nc <- nc_open(fname)
  
  x <- ncvar_get(nc, var)
  idx <- apply(
    x,
    2,
    function(x) sum(
      !is.na(x)
    )
  )
  idx <- idx != 0
  
  obs <- x[ ,idx]
  
  nS <- nrow(Xf)
  nM <- ncol(Xf)
  
  Xd <- matrix(list(), nS, nM)
  
  for (i in 1:nS) {
    for (j in 1:nM) {
      
      fut <- Xf[[i, j]]
      cont <- Xc[[j]]
      
      idx <- !(is.null(obs) | is.null(cont) | is.null(fut))
      
      if (idx) {
        cont <- cont*mult
        fut <- fut*mult
        x <- Downscaling(obs, cont, fut, A2)
        Xd[[i, j]] <- x
      }
      
    }
  }
  return(Xd)
}

Downscaling <- function(observation, control, future, var) {
  
  # input
  
  #      observation: vector/matrix of observations(if different stations rows:days x columnS:Nstation) same length of control
  #      control: vector/matrix of RCM output for control(if different stations rows:days x columnS:Nstation) same length of obs
  #      future: vector/matrix of RCM output for future(if different stations rows:days x columnS:Nstation)
  #      var: 'P' for precipitation,
  #           'T' for temperature
  
  
  # output
  
  #           ControlD: vector/matrix of the RCM output downscaled
  #           FutureD: vector/matrix of RCM future output downscaled
  
  
  # settings mto be done here:
  #                              f : semi-amplitude of the moving window
  #                              TP: precipitation threshold
  
  
  # functions:
  #           compute_downscaling_qq_Daily_MW
  #           apply_downscaling_qq_Daily_MW
  
  
  
  switch(var,
         "T"={A1 = ' Temperature [C?] ';
         }, "P"={A1 = ' Precipitation [mm/d] ';},
         stop("Variable name is undefined"))
  
  Ty= 365; 
  
  ControlD <- matrix(NA, nrow(control) - Ty*2, ncol(control));
  FutureD = matrix(NA, nrow(future), ncol(future))
  for (ss in 1:ncol(observation)){ #--> ss= #stations
    
    y_obs = observation[ ,ss] ;
    y_ctrl = control[ ,ss] ;
    
    
    #plot( y_obs, col = "red", xlab='Time [days]', ylab = A1, type = "l")
    #lines( y_ctrl, col = 'blue' )
    
    # CFD
    obs_cum_fun <-  ecdf(y_obs); ycdf_obs <- obs_cum_fun(y_obs);
    #plot(y_obs, ycdf_obs, col = 'red', ylab = 'Empirical cdf', xlab = A1)
    ctrl_cum_fun <-  ecdf(y_ctrl); ycdf_ctrl <- ctrl_cum_fun(y_ctrl);
    #points(y_ctrl, ycdf_ctrl, col = 'blue')
    
    ## SETTINGS
    w=45; # days # semi-aplitude of the moving window
    TP=1; # precipitation treshold mm/d
    
    ## DOWNSCALING - calibration phase
    
    y_ctrl_d <- compute_downscaling_qq_Daily_MW(y_obs , y_ctrl , var , w , 1 ,TP) ;
    down.param <- y_ctrl_d[[1]]; y_ctrl_d <- y_ctrl_d[[2]] 
    #[ y_ctrl_d, dwnsc_param , H ] = compute_downscaling_qq_Daily_MW(y_obs , y_ctrl , var , w , 1 ,TP) ;
    
    ControlD[ ,ss] <- y_ctrl_d;
    
    # need to remove the first and last year  used to build the moving window
    # for checking
    
    y_obs_short <- y_obs[Ty:(length(y_obs)-Ty-1)]; y_ctrl_short <- y_ctrl[Ty:(length(y_ctrl)-Ty-1)];
    
    # check1: ctrl trajectory after downscaling
    
    #plot(y_obs_short, col = 'black', xlab = ('Days'), ylab = (A1))
    #lines( y_ctrl_short, col = 'red' )
    #lines( y_ctrl_d, col = 'cyan' )
    
    # # check2: if you downscale y_bck --> qq plot must be a line 45 degrees
    # 
    pc = quantile(y_ctrl_d[y_ctrl_d>TP] , 1:99/100 ) ;
    po = quantile( y_obs_short[y_obs_short>TP] , 1:99/100 ) ;
    Lim = max( max(po) , max(pc) ) ;
    #plot( pc,po, xlim = c(0, Lim), ylim  = c(0, Lim), type = "l")
    #lines(c(0, Lim), c(0, Lim), col = "red")
    
    ## DOWNSCALING - projection phase
    y_scen = future[ ,ss];
    
    y_scen_d  = apply_downscaling_qq_Daily_MW(y_scen ,var, down.param[[4]] , down.param[[3]], TP);
    
    # check: scen trajectory after downscaling
    #plot( y_scen, col = 'magenta', xlab = ('Days'), ylab = (A1))
    #lines( t(y_scen_d), col = 'cyan' )
    FutureD[ ,ss] = t(y_scen_d) ;
  }
  return(list(ControlD, FutureD))
}    

compute_downscaling_qq_Daily_MW <- function( obsr , ctrl ,  var , w ,  picture, TP ) {
  
  # [ ctrl_d , dwnsc_param , H ] = compute_downscaling_qq_Daily_MW( obsr , ctrl ,  var , f , picture, TP)
  #
  # The function downscales precipitation or temperature series simulated
  # through a climate model with respect to observed data.
  # Downscaling is due to the quantile-quantile (Q-Q) plot between 
  # observed time-serie and control time-serie, simulated through a climate
  # model over the same time period of the observed serie.
  # 
  # A Q-Q plot is a probability plot, a kind of graphical method for 
  # comparing two probability distributions, by plotting their quantiles
  # against each other. The Q-Q plot is used as a correction function
  # in order to match the ctrl cdf (pdf) with the the obrs cdf (pdf).
  # The scen time-serie, i.e. the serie simulated through the climate model
  # for a future period, is downscaled using of the same correction function 
  # used to downscale the ctrl time serie: we ignore the cdf of the future
  # climate, and the most reasonable hypothesis is that the climate model's
  # errors will be the same as in the past day climate.
  # 
  # The correction function is based on the 99 percentiles of the ctrl and 
  # obsr cdfs. Outside this range (i.e. < 1 prc or > 99 prc), a constant 
  # correction is applied: if the 99th quantile of January temperature is 
  # corrected by +1?C, any temperature above this threshold is corrected by 
  # +1?C. A linear interpolation is applied between two percentiles.
  # 
  # The correction function is estimated on a daily (365 correction 
  # functions only) with a moving window.
  #
  # Input:
  #   - obsr  = observed time-serie                                - vector [ h , 1 ]
  #   - ctrl  = control time-serie                                 - vector [ h , 1 ]
  #   - var         = type of variable to be downscaled.           - string
  #                   It can assume two values only: 'P' = precipitation
  #                                                  'T' = temperature
  #   - w     = semi-amplitude of the window - scalar                                                 
  #   - picture     = flag (0,1)                                   - scalar   
  #                   If different from zero , Q-Q plot is created
  # Output:
  #   - ctrl_d      = downscaled control time-serie                - vector [ h_ctrl , 1 ]
  #   - dwnsc_param = tables with the percentiles of the 
  #                   observation and control time series          - struct
  #   - H           = handles to lineseries graphics objects       - vector [ I , 1]
  
  # time period
  dwnsc_param <- matrix(list(), 4, 1)
  Ty=365;  # days
  
  # According to the semi-amplitude of the window (f), the
  # first and last years of the time series are used only to build the first and
  # last windows for the first/last f days. f corresponds to the
  # semi-amplitude of the windows to center the windows on the first day of a window with amplitude 2 x f. 
  # Observation and control vectors have to be reshaped in order to remove
  # the first/last year and add the f days at the beginning and at the end
  # of the time series. 
  
  N_years  = length( obsr ) / Ty;
  Y_obs = obsr; dim(Y_obs) <- c(Ty, N_years)
  Y_ctrl= ctrl; dim(Y_ctrl) <- c(Ty, N_years);
  
  # extraction of the last f days of December of the first year 
  y_ad_top_oss=Y_obs[(nrow(Y_obs)-w+1):nrow(Y_obs), 1];
  y_ad_top_ctrl=Y_ctrl[(nrow(Y_obs)-w+1):nrow(Y_obs), 1];
  
  # extraction of the first f days of January of the last year 
  y_ad_down_oss=Y_obs[1:w, ncol(Y_obs)]; y_ad_down_oss <- y_ad_down_oss[length(y_ad_down_oss):1]
  y_ad_down_ctrl=Y_ctrl[1:w, ncol(Y_ctrl)]; y_ad_down_ctrl <- y_ad_down_ctrl[length(y_ad_down_ctrl):1]
  
  # removal of the first and the last year
  obsr=Y_obs[ ,2:(ncol(Y_obs)-1)];
  dim(obsr) <- c(nrow(obsr)*ncol(obsr), 1)
  y_obs_original=obsr;
  
  ctrl=Y_ctrl[ ,2:(ncol(Y_ctrl)-1)];
  dim(ctrl) <- c(nrow(ctrl)*ncol(ctrl), 1)
  y_ctrl_original=ctrl;
  
  # Addition of the first and last f days
  obsr=c(y_ad_top_oss, obsr , y_ad_down_oss);
  ctrl=c(y_ad_top_ctrl, ctrl, y_ad_down_ctrl);
  
  switch(var,
         "T"={C_I = -50; C_S = 50;
         }, "P"={C_I = 0; C_S = 500;},
         stop("Variable name is undefined"))
  
  
  # matrix  [ h , 365 ] of logical value: 1s in the j-th column stand for
  # data of the time serie belonging to the j-th day of the year. The matrix
  # ID is the index of the day of the year on which the correction function
  # is estimated. 
  
  Ny = floor( length(y_obs_original)/365 ) ;
  ID <- matrix( rep( diag(1, 365) , Ny ) , ncol =  ncol(diag(1, 365)) , byrow = TRUE )
  I=ncol(ID); 
  
  # percentiles 
  Prc = (1:99)/100 ;
  
  # Prc_obsr and Prc_ctrl: two correction table with the 99 percentiles
  # of the two distributions (obsr data and simulation on the ctrl period)        
  
  Prc_obsr <- matrix(NA, length(Prc), I)
  Prc_ctrl <- matrix(NA, length(Prc), I)
  for (i in 1 : I){
    
    # For each day of the year i, the correction function is estimated
    # using 2 x f + 1 values centered in the day i es?xtracting the f days before and f days after i. 
    
    
    # The same day i are  found in the time series.   
    
    k <- ID[ ,i];  k = which(k == 1); #[# years,1]  
    
    # for example if my time series is 50 years long, I'm going to find 50 1-Jan in k.    
    # Around each 1-Jan I need to build the window: [f days - 1 Jan - f days]
    
    ob <- list()
    ctr <- list()
    for (z in 1:length(k)){
      ob[[z]] <- obsr[k[z]:(k[z]+2*w)];
      ctr[[z]]=ctrl[k[z]:(k[z]+2*w)];
    }
    ob <- unlist(ob); ctr <- unlist(ctr)
    # Thus, creating the values for the percentile estime for the i day
    
    switch(var,
           "T"={Prc_temp = as.numeric(quantile( ob , Prc )) ; 
           Prc_obsr[ ,i] = t(Prc_temp) ;
           Prc_ctrl_temp = as.numeric(quantile( ctr , Prc )) ; 
           Prc_ctrl[ ,i] = t(Prc_ctrl_temp);
           }, "P"={ob[ob < TP ] = 0;
           Prc_temp = as.numeric(quantile( ob[ob>TP] , Prc )) ; 
           Prc_obsr[ ,i] = t(Prc_temp) ;
           Prc_ctrl_temp = as.numeric(quantile( ctr[ctr>TP] , Prc )) ; 
           Prc_ctrl[ ,i] = t(Prc_ctrl_temp) ;},
           stop("Variable name is undefined"))
  }
  
  # # FIGURE: QQ plot
  # if picture
  # figure
  # switch var
  # case 'P'
  # Lim = max( max(Prc_obsr(end,:)) , max(Prc_ctrl(end,:)) ) ;
  # H = plot( Prc_ctrl , Prc_obsr ) ;
  # #             opz_line = { '-' , ':' , '--' , '-.' } ;
  # #             for i = 1:I
  # #                 set( H(i) , 'LineWidth' , 1 , 'LineStyle' , opz_line{i} )
  # #             end
  # hold on; plot( [0 Lim] , [0 Lim] , 'Color' , [.5 .5 .5] )
  # set( gca , 'XLim' , [0  Lim] , 'YLim' , [0  Lim] )
  # 
  # xlabel( 'Backcast' )
  # ylabel( 'Observation' )
  # grid on
  # case 'T'
  # Lim1 = min( min(Prc_obsr(1,:)) , min(Prc_ctrl(1,:)) ) ;
  # Lim2 = max( max(Prc_obsr(end,:)) , max(Prc_ctrl(end,:)) ) ;
  # H = plot( Prc_ctrl , Prc_obsr ) ;
  # #             opz_line = { '-' , ':' , '--' , '-.' } ;
  # #             for i = 1:I
  # #                 set( H(i) , 'LineWidth' , 1 , 'LineStyle' , opz_line{i} )
  # #             end
  # hold on; plot( [Lim1 Lim2] , [Lim1 Lim2] , 'k' )
  # set( gca , 'XLim' , [Lim1  Lim2] , 'YLim' , [Lim1  Lim2] )
  # xlabel( 'Backcast' )
  # ylabel( 'Observation' )
  # grid on
  # end
  # else
  #   H = [] ;
  # end
  
  dwnsc_param[[1]] = Prc_ctrl ;
  dwnsc_param[[2]] = Prc_obsr ;
  # # -----------------------------------------------------------------------
  # # figure useful to understand the meaning of the following variables:
  # Prc_C = sort( rand(99,1) -0.5 ) ;
  # Prc_O = sort( rand(99,1) - 0.5) ;
  # figure; plot( Prc_C , Prc_O , '.' )                                         # Q-Q plot
  # hold on 
  # plot( min(Prc_C) , min(Prc_O) , 'oc' )                                      # min
  # plot( max(Prc_C) , max(Prc_O) , 'ok' )                                      # max
  # plot( -1 , -1 + ( min(Prc_O) - min(Prc_C) ) , 'sc')                         # inf
  # plot( 1 , 1 + ( max(Prc_O) - max(Prc_C) ) , 'sk')                           # sup
  # plot( [-1 min(Prc_C)],[ -1 + ( min(Prc_O) - min(Prc_C) ) min(Prc_O)], '-c') # constant correction
  # plot( [1 max(Prc_C)] , [1 + ( max(Prc_O) - max(Prc_C) ) max(Prc_O)] , '-k') # constant correction
  # plot( [-1 1] , [-1 1] , 'r' )     
  # legend( 'Q-Q plot' , 'min' , 'max' , 'inf' , 'sup' )
  # xlabel('ctrl'); ylabel('obsr')
  # # -----------------------------------------------------------------------
  
  # Outside the table range (< 1 prc or > 99 prc), a constant correction is
  # applied:
  # e.g. if the 99th quantile of January temperature is corrected by +1?C,
  # any temperature above this threshold is corrected by +1?C     
  ctrl_min = apply(Prc_ctrl, 2, min)                       ;
  ctrl_max = apply(Prc_ctrl,2,max)                     ;
  ctrl_inf <- matrix( rep( C_I , ncol(Prc_ctrl) )) ;
  ctrl_sup <- matrix( rep( C_S , ncol(Prc_ctrl) )) ;
  Prc_ctrl = rbind(t(ctrl_inf), Prc_ctrl, t(ctrl_sup) )    ;  # [ length(Prc)+2 , I ]
  
  obsr_min = apply(Prc_obsr, 2 ,min)                    ;
  obsr_max = apply(Prc_obsr, 2, max)                    ;
  obsr_inf = ctrl_inf + ( obsr_min - ctrl_min ) ;     
  
  
  switch(var,
         "P"={obsr_inf[obsr_inf <0] = 0 ;})
  
  obsr_sup = ctrl_sup + ( obsr_max - ctrl_max ) ; 
  Prc_obsr = rbind(t(obsr_inf), Prc_obsr, t(obsr_sup) ) ;     # [ length(Prc)+2 , I ]
  
  
  # A linear interpolation is applied between two percentiles: 
  # the values of Prc_ctrl should be distinct, if they are not 
  # (e.g. dry climate) an eps will be added.
  
  
  ctrl_d = matrix(0, nrow(y_ctrl_original), ncol(y_ctrl_original));
  
  for (i in 1:I){
    
    
    #Prc_ctrl( find( diff( Prc_ctrl(:,i) ) == 0  )+1 , i ) = Prc_ctrl( find( diff( Prc_ctrl(:,i) ) == 0  ) , i ) + eps .* find( diff( Prc_ctrl(:,i) ) == 0 ) ;
    # data to be modifies (x' = previous)        = data before the one to be modified                 + eps .* 99 (max, in case of 98 percentile  equal to 99 percentile, very difficult)
    
    #[~, uidx] = unique(Prc_ctrl(:,i)); #  The grid vectors must contain unique points
    library(pracma)
    switch(var,
           "T"={idx_temp1=as.logical(ID[,i]);
           }, "P"={temp1=y_ctrl_original*ID[ ,i];
           temp1[temp1<TP]=0;
           temp1[temp1>=TP]=1;
           idx_temp1=as.logical(temp1);},
           stop("Variable name is undefined"))
    inp_mod <- data.frame(x =  Prc_ctrl[ ,i], y = Prc_obsr[ ,i]) 
    mod <- lm(y~x, data = inp_mod)
    ctrl_d[idx_temp1]= predict(mod, newdata= data.frame(x = y_ctrl_original[idx_temp1]))
  }
  
  # Vectors of downscaled time-series [ h , 1 ]:
  dwnsc_param[[3]] = Prc_ctrl ;
  dwnsc_param[[4]] = Prc_obsr ;
  return(list(dwnsc_param, ctrl_d))
}

apply_downscaling_qq_Daily_MW <- function( scen ,var, Prc_obsr , Prc_ctrl, TP) {
  
  
  # [ scen_d ] = apply_downscaling_qq( scen, tdates_scen, Prc_obsr, Prc_ctrl)
  #
  # The function downscales precipitation or temperature timeseries simulated
  # through a climate model. Downscaling is made by the Quantile Method. 
  # 
  # The downscaling function is given by a Q-Q plot estimated on a daily (365 
  # correction functions) computed on a moving window. The window has an
  # amplitude 2 x f and is centered in each day of the year.
  # A Q-Q plot is a probability plot, a kind of graphical method for 
  # comparing two probability distributions, by plotting their quantiles
  # against each other. 
  # The correction function is based on the 99 percentiles of the 
  # distribution of the observation and of the control timeserie.
  # A linear interpolation is applied between two percentiles.
  # Outside that range (i.e. < 1 prc or > 99 prc), a constant 
  # correction is applied: if the 99th quantile of temperature is corrected
  # by +1?C, any temperature above this threshold is corrected by +1?C. 
  # 
  #
  # Input:
  #   - scen        = scenario time-serie                          - vector [ h_scen , 1 ] 
  #   - var         = type of variable to be downscaled.           - string
  #                   It can assume two values only: 'P' = precipitation
  #                                                  'T' = temperature
  #   - Prc_obsr    = percentiles of the observation               - matrix [101,I]
  #   - Prc_ctrl    = percentiles of the control timeserie         - matrix [101,I]
  #                   
  # Output:
  #   - scen_d      = downscaled scenario timeserie                - vector [ h_scen , 1 ]
  # 
  
  
  
  
  # matrix  [ h , 365 ] of logical value: 1s in the j-th column stand for
  # data of the time serie belonging to the j-th day of the year. The matrix
  # ID is the index of the day of the year for which the correction function
  # is estimated. 
  
  I = ncol(Prc_ctrl);
  Ny = floor(length(scen)/365 ) ;
  ID = matrix( rep( diag(1, 365) , Ny ) , ncol =  ncol(diag(1, 365)) , byrow = TRUE )
  scen_d = scen*0;
  
  # A linear interpolation is applied between two percentiles: 
  # the values of Prc_ctrl should be distinct, if they are not 
  # (e.g. dry climate) an eps will be added.
  for (i in 1:I){
    uidx = as.numeric(rownames(unique(data.frame(Prc_ctrl[ ,i])))) 
    
    switch(var,
           "T"={idx_temp1 <- as.logical(ID[ ,i]);
           }, "P"={temp1=scen*ID[ ,i] ;
           temp1[temp1<TP]=0;
           temp1[temp1>=TP]=1;
           idx_temp1=as.logical(temp1);
           warnings()},
           stop("Variable name is undefined"))
    
    #     scen_d( find(ID(:,i)) ) = interp1( Prc_ctrl(uidx,i) , Prc_obsr(uidx,i) , scen( find(ID(:,i)) ) ) ;
    inp_mod <- data.frame(x =  Prc_ctrl[uidx ,i], y = Prc_obsr[uidx ,i]) 
    mod <- lm(y~x, data = inp_mod)
    scen_d[idx_temp1]= predict(mod, newdata= data.frame(x = scen[idx_temp1])); 
  }
  
  # Vectors of downscaled time-series [ h , 1 ]:
  scen_d = t(scen_d) ;
  scen_d[scen_d < 0] <- 0
  return(scen_d)
}

nc_from_ts <- function(x, init, var, id, pss) {
  
  # nc_from_ts <- function(x, init, var, id, pss) 
  # convert time series to netcdf files
  # Input:
  #
  # x = time series file
  # init = path to folder 
  # var = variable name
  # id = model
  # pss = scenario
  #
  # Returns:
  # NA
  
  fin <- paste(init, var, sep = '/')
  setwd(fin)
  
  fname <- list.files(pattern = '*.nc')
  name <- grepl('obs', fname)
  
  
  nc <- nc_open(fname[name])
  
  lon <- ncvar_get(nc, 'lon')
  lat <- ncvar_get(nc, 'lat')
  
  y <- ncvar_get(nc, var)
  val <- apply(
    y,
    2,
    function(x) sum(
      !is.na(x)
    )
  )
  nc_close(nc)
  val <- val != 0
  
  nR <- nrow(x)
  nC <- ncol(x)
  for (r in 1:nR) {
    for (c in 1:nC) {
      
      cont <- x[[r, c]][[1]]
      fut <- x[[r, c]][[2]]
      
      idx <- !(is.null(cont) | is.null(fut))
      
      if (idx) {
        
        Dc <- matrix(NA, nrow(cont), ncol(y))
        Df <- matrix(NA, nrow(fut), ncol(y))
        
        Dc[ ,val] <- cont
        Df[ ,val] <- fut
        
        nc_gen(Dc, var, 'control', id[c], lat, lon)
        nc_gen(Df, var, pss[r], id[c], lat, lon)
      }
      
    }
  }
}

nc_gen <- function(z, var, label_x, label_y, lat, lon) {
  
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
  
  dim_idx <- ncdim_def(name='Index',
                       units='m',
                       longname='idx',
                       vals= 1:ncol(z) )
  
  dim_time <- ncdim_def('time',
                        units='days from 2005-12-31',
                        longname='time',
                        calendar="standard", vals=1:nrow(z)
  )
  
  varLat <- ncvar_def(name='lat',
                      units='degrees_north',
                      dim=list(dim_idx),
                      missval=NA,
                      longname='latitude',
                      prec='double'
  )
  
  varLon <- ncvar_def(name='lon',
                      units='degrees_east',
                      dim=list(dim_idx),
                      missval=NA,
                      longname='latitude',
                      prec='double'
  )
  
  varX <- ncvar_def(name=A,
                    units= u,
                    dim=list(dim_time, dim_idx),
                    missval=NA,
                    longname=A1
  )
  
  vars <- list(varLat, varLon, varX)
  
  outputfile <- paste('Umbeluzi', label_x, label_y, '.nc', sep = '_')
  
  con <- nc_create(outputfile, vars)
  ncvar_put(con, varLat, lat)
  ncvar_put(con, varLon, lon)
  ncvar_put(con, varX, z)
  
  nc_close(con)
  
}

# Postprocess----

futures_plot_data <- function(init_n, var, labs, id) {
  
  nS <- length(labs)
  nM <- length(id)
  rm <- as.numeric(
    as.Date('2040-01-01')
    - 
      as.Date('2006-01-01')
  )
  
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
  
  x <- matrix(NA, nS, nM)
  
  for (i in 1:nS) {
    f2 <- paste(var, labs[i], sep = '/')
    setwd(paste(init_n, f2, sep = '/'))
    
    f <- list.files(pattern = '*.nc')
    for (j in 1:nM) {
      idx <- grepl(id[j], f)
      
      if (sum(idx) != 0) {
        nc <- nc_open(f[idx])
        v <- ncvar_get(nc, A)
        if (i != 1) {
          v <- v[-(1:rm),]
        }
        x[i, j] <- mean(v, na.rm = T)
      }
      
    }  
  }
  return(x)
}

plot_futures <- function(p, vars, pss, id) {
  
  i <- 1
  
  for (var in vars) {
    pi <- p[[i]]
    
    y <- matrix(NA, nrow(pi)-1, ncol(pi))
    
    for (r in 1:nrow(y)) {
      y[r, ] <- pi[r+1, ] - pi[1, ]
    }
    
    y <- as.data.frame(y)
    y <- cbind(y, c('4.5', '8.5'))
    colnames(y) <- c(id, 'RCP')
    
    y <- melt(y, 'RCP')
    
    if(i == 1) {
      z <- y
    }else{
      z <- cbind(z, y[ ,ncol(y)])
    }
    i <- i + 1
  }
  colnames(z)[-1] <- c('Model', vars)
  x0 <- seq(
    min(
      min(
        z[ ,3],
        na.rm  = T),
      0
    ),
    max(
      z[ ,3],
      na.rm = T
    ),
    length.out =
      nrow(z)
  )
  
  y0 <- seq(min(
    min(
      z[ ,4],
      na.rm  = T),
    0
  ),
  max(
    z[ ,4],
    na.rm = T),
  length.out =
    nrow(z)
  )
  
  z0 <- rep(0, nrow(z))
  b <- data.frame(x0, y0, z0)
  
  ggplot()+
    geom_point(data = z,
               aes(pav,
                   tas,
                   col = RCP,
                   shape = Model),
               size = 4)+
    xlab(
      bquote(
        Delta ~ 'P' ~ '[mm/day]')
    )+
    ylab(
      bquote(Delta ~ 'T' ~ '[K]')
    )+
    scale_color_manual(
      values = c('green',
                 'red',
                 'blue')
    )+
    ggtitle('Anomalies 2040-2100')+
    geom_line(data = b,
              aes(z0, y0),
              col = 'black')+
    geom_line(data = b, aes(x0, z0), col = 'black')+
    theme_bw()#+
    #theme(text=element_text(size=12,  family="Calibri"))
  
  ggsave('Anomalies.pdf',
         plot = last_plot(),
         device = 'pdf',
         scale = 1, 
         dpi = 300)
}

MASH_data <- function(var, init_n, model, periods) {
  
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

MASH_plot <- function(x, Ys, w, Y, datalab, titolo) {
  #
  # function [m,h]= MASH_plot(x, Ys, w, Y, datalab) 
  #
  # Computes and displays the Moving Average over Shifting Horizon (MASH), 
  # a tool to be applied in Exploratory Data Analysis for trend detection
  # in seasonal data. 
  # The MASH allows for simultaneously investigating the seasonality in the 
  # data and filtering out the effects of interannual variability,
  # thus making trend detection easier.
  # Ref. "Trend detection in seasonal data: from hydrology to water
  # resources", D. Anghileri, F. Pianosi, and R. Soncini-Sessa, Journal of 
  # Hydrology, 2014 (http://dx.doi.org/10.1016/j.jhydrol.2014.01.022)
  #
  # The MASH is a collection of trajectories of the moving average over 'w' 
  # consecutive days (or weeks, months, etc.), computed on a shifting horizon
  # of 'Y' consecutive years.
  #
  # Input:
  #    - x  = time series to be analyzed (*)                    - matrix(T,nY)
  #    - Ys = starting year of the time series (format: 'yyyy') - scalar
  #    - w  = half number of consecutive data to be averaged    - scalar
  #          (i.e. moving average considers a window of 2*w+1 data)
  #    - Y  = number of consecutive years to be averaged (**)   - scalar
  # datalab = name and units of measure of 'x' (***)            - string
  #
  # Output:
  #   - m  = averaged seasonal data                            - matrix(T,nY-Y+1)
  #   - h  = handle of the lines in the plot                   - vector(nY-Y+1,1)
  #
  # REMARKS:
  # (*) 'T' is the number of calendar units in one period (e.g. 365 for daily
  #      data, 52 for weekly data, etc.) and 'nY' is the number of years in 
  #      the time series. In case of daily data, data of leap years must be
  #      resorted to a time series of 365 values by the user.
  # (**) If Y = 1, w = 0 then the original data are plotted (without
  #      averaging)
  # (***) To be used for labeling axes in the plots
  #
  #
  # Created:     Daniela Anghileri, Francesca Pianosi 28/01/2014
  
  library(reshape2)
  library(ggplot2)
  
  # check on inputs
  if (!is.numeric(x)) { stop('x must be a double matrix of size (Tp,nY)')}
  Tp <- dim(x)[1]; nY <- dim(x)[2] ;
  if (!is.numeric(Ys)) { stop('Ys must be a scalar integer')}
  if (!is.atomic(Ys)) { stop('Ys must be a scalar integer')}
  if (!is.numeric(Y)) { stop('Y must be a scalar integer')}
  if (!is.atomic(Y)) { stop('Ys must be a scalar integer')}
  
  
  if (Y>nY){stop(paste('Y must be <= '),nY)} else if (Y<1) {stop('Y must be >=1')}
  if (!is.numeric(w)) { stop('w must be a scalar integer')}
  if (!is.atomic(w)) { stop('w must be a scalar integer')}
  
  if (2*w+1 > Tp){
    stop(paste('w must be such that 2*w+1 be <= #d (number of rows in x)',Tp) )
  }else if (w<0){stop('w must be >=0')}
  
  if (!is.character(datalab)){stop('datalab must be a string')}
  
  # define useful variables
  Hi   = c(Ys:(Ys+nY-Y))   ;
  H_if = c(Hi, Hi+(Y-1)) ; # initial and final year
  
  
  # Compute the MASH:
  m  = matrix(NA, Tp, nY-Y+1);
  x_ = rbind(x[(nrow(x)-w+1):nrow(x), ], x, x[1:w, ] ) ; # handling first and last w data
  for (j in 1: (nY-Y+1) ){
    for (i in 1:Tp){
      m[i,j] = mean( x_[ i:(i+2*w) , j:(j+Y-1) ], na.rm = T)
    }
  }
  m <- t(m); mash <- m
  
  
  # Prepare data for plot
  m <- data.frame(Year = Hi, m)
  colnames(m)[2:ncol(m)] <- c(1:365)
  m <- melt(m, "Year"); m[ ,2] <- as.numeric(m[ ,2])
  
  Months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sept", "Oct", "Nov", "Dec")
  
  m_p <- m[m[ ,1] <= 2005, ]
  
  m_min <- m_p %>% group_by(variable) %>%
    slice(min(value))
  
  m_min <- aggregate(value~variable, m_p, min)
  m_max <- aggregate(value~variable, m_p, max)
  m_s <- cbind(m_min, m_max[ ,2])
  colnames(m_s)[3] <- 'v2' 
  
  # Plot
  ggplot(data = m_p, aes(variable, value, group = Year,  col  = Year))+
    geom_line()+
    scale_color_gradient2(low = "springgreen",
                          mid = "yellow",
                          high = "red3", 
                          midpoint = (max(m_p[ ,1]) + min(m_p[ ,1])) /2 ,
                          limits = c(min(m_p[ ,1]), max(m_p[ ,1]) ) )+
    xlab("Month")+
    ylab(datalab)+
    ggtitle(paste('Hystorical', titolo))+
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5)#,
          #text=element_text(family="Calibri")
          )+
    scale_x_continuous(breaks = round(seq(min(m[ ,2]), max(m[ ,2]), by = 365/12), 0), 
                       labels = Months)
  
  ggsave(paste0("H_Mash",var, ".pdf"),
         plot = last_plot(), device = "pdf", scale = 1, dpi = 300)
  
  #return(mash)
  
  m_p <- m[m[ ,1] >= 2005, ]
  
  ggplot()+
    geom_ribbon(data = m_s, aes(x = variable, ymin = value, ymax = v2), alpha = 0.8 )+
    geom_line(data = m_p, aes(variable, value, group = Year,  col  = Year))+
    scale_color_gradient2(low = "springgreen",
                          mid = "yellow",
                          high = "red3", 
                          midpoint = (max(m_p[ ,1]) + min(m_p[ ,1])) /2 ,
                          limits = c(min(m_p[ ,1]), max(m_p[ ,1]) ) )+
    xlab("Month")+
    ylab(datalab)+
    ggtitle(titolo)+
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5)#,
          #text=element_text(family="Calibri")
          )+
    scale_x_continuous(breaks = round(seq(min(m[ ,2]), max(m[ ,2]), by = 365/12), 0), 
                       labels = Months)
  
  ggsave(paste0("Mash",var, ".pdf"),
         plot = last_plot(), device = "pdf", scale = 1, dpi = 300)
  
  
  
} 



