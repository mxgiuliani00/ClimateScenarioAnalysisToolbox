function [ ctrl_d , dwnsc_param , H ] = compute_downscaling_qq_Daily_MW( obsr , ctrl ,  var , w ,  picture, TP )

% [ ctrl_d , dwnsc_param , H ] = compute_downscaling_qq_Daily_MW( obsr , ctrl ,  var , f , picture, TP)
%
% The function downscales precipitation or temperature series simulated
% through a climate model with respect to observed data.
% Downscaling is due to the quantile-quantile (Q-Q) plot between 
% observed time-serie and control time-serie, simulated through a climate
% model over the same time period of the observed serie.
% 
% A Q-Q plot is a probability plot, a kind of graphical method for 
% comparing two probability distributions, by plotting their quantiles
% against each other. The Q-Q plot is used as a correction function
% in order to match the ctrl cdf (pdf) with the the obrs cdf (pdf).
% The scen time-serie, i.e. the serie simulated through the climate model
% for a future period, is downscaled using of the same correction function 
% used to downscale the ctrl time serie: we ignore the cdf of the future
% climate, and the most reasonable hypothesis is that the climate model's
% errors will be the same as in the past day climate.
% 
% The correction function is based on the 99 percentiles of the ctrl and 
% obsr cdfs. Outside this range (i.e. < 1 prc or > 99 prc), a constant 
% correction is applied: if the 99th quantile of January temperature is 
% corrected by +1°C, any temperature above this threshold is corrected by 
% +1°C. A linear interpolation is applied between two percentiles.
% 
% The correction function is estimated on a daily (365 correction 
% functions only) with a moving window.
%
% Input:
%   - obsr  = observed time-serie                                - vector [ h , 1 ]
%   - ctrl  = control time-serie                                 - vector [ h , 1 ]
%   - var         = type of variable to be downscaled.           - string
%                   It can assume two values only: 'P' = precipitation
%                                                  'T' = temperature
%   - w     = semi-amplitude of the window - scalar                                                 
%   - picture     = flag (0,1)                                   - scalar   
%                   If different from zero , Q-Q plot is created
% Output:
%   - ctrl_d      = downscaled control time-serie                - vector [ h_ctrl , 1 ]
%   - dwnsc_param = tables with the percentiles of the 
%                   observation and control time series          - struct
%   - H           = handles to lineseries graphics objects       - vector [ I , 1]

% time period

T=365;  % days

% According to the semi-amplitude of the window (f), the
% first and last years of the time series are used only to build the first and
% last windows for the first/last f days. f corresponds to the
% semi-amplitude of the windows to center the windows on the first day of a window with amplitude 2 x f. 
% Observation and control vectors have to be reshaped in order to remove
% the first/last year and add the f days at the beginning and at the end
% of the time series. 

y_obs_original=obsr;
y_ctrl_original=ctrl;
N_years  = length( obsr ) / T  ;
Y_obs= reshape( obsr , T , N_years ) ;
Y_ctrl= reshape(ctrl , T , N_years ) ;

% extraction of the last f days of December of the first year 
y_ad_top_oss=Y_obs(end-w+1:end, 1);
y_ad_top_ctrl=Y_ctrl(end-w+1:end, 1);

% extraction of the first f days of January of the last year 
y_ad_down_oss=flip(Y_obs(1:w, end));
y_ad_down_ctrl=flip(Y_ctrl(1:w, end));

% Addition of the first and last f days
obsr=[y_ad_top_oss;obsr  ;y_ad_down_oss];
ctrl=[y_ad_top_ctrl;ctrl ;y_ad_down_ctrl];


switch var
    
    case 'P'
        C_I =   0 ; % mm/day
        C_S = 500 ; % mm/month--> 300/30 mm/day
        
    case 'T'
        C_I = -50 ; % °C
        C_S =  50 ; % °C
       
end

% matrix  [ h , 365 ] of logical value: 1s in the j-th column stand for
% data of the time serie belonging to the j-th day of the year. The matrix
% ID is the index of the day of the year on which the correction function
% is estimated. 

Ny = floor( length(y_obs_original)/365 ) ;
ID = repmat( eye(365), Ny, 1 ) ;
I=size(ID,2); 


% percentiles 
Prc = 1:99 ;

% Prc_obsr and Prc_ctrl: two correction table with the 99 percentiles
% of the two distributions (obsr data and simulation on the ctrl period)        
     
  
for i = 1 : I
   
    ob=[];
    ctr=[];
    
% For each day of the year i, the correction function is estimated
% using 2 x f + 1 values centered in the day i esìxtracting the f days before and f days after i. 

    
% The same day i are  found in the time series.   

    k=find(ID(:,i)); %[# years,1] 
    
% for example if my time series is 50 years long, I'm going to find 50 1-Jan in k.    
% Around each 1-Jan I need to build the window: [f days - 1 Jan - f days]
    
    
         for z=1:length(k)
        
         ob_temp=obsr( k(z):k(z)+2*w);
         ob=[ob;ob_temp]; 
         ctr_temp=ctrl( k(z):k(z)+2*w);
         ctr=[ctr;ctr_temp];  
         
         end

% Thus, creating the values for the percentile estime for the i day
         

 
   switch var
       
       case 'P'
           
 % In case of P, data for dirzzling need to be corrected. Precipitation
 % below 1 mm/day is set to zero and the rest above the threshold is used
 % to estime percentile. 
 
    ob(ob < TP ) = 0;
    Prc_temp = prctile( ob(ob>TP) , Prc ) ; 
    Prc_obsr(:,i) = Prc_temp' ;
    Prc_ctrl_temp = prctile( ctr(ctr>TP) , Prc ) ; 
    Prc_ctrl(:,i) = Prc_ctrl_temp' ;
    
       case 'T'
           
    
    Prc_temp = prctile( ob , Prc ) ; 
    Prc_obsr(:,i) = Prc_temp' ;
    Prc_ctrl_temp = prctile( ctr , Prc ) ; 
    Prc_ctrl(:,i) = Prc_ctrl_temp' ;
    
    end
   
end

% FIGURE: QQ plot
if picture
    figure
    switch var
        case 'P'
            Lim = max( max(Prc_obsr(end,:)) , max(Prc_ctrl(end,:)) ) ;
            H = plot( Prc_ctrl , Prc_obsr ) ;
%             opz_line = { '-' , ':' , '--' , '-.' } ;
%             for i = 1:I
%                 set( H(i) , 'LineWidth' , 1 , 'LineStyle' , opz_line{i} )
%             end
            hold on; plot( [0 Lim] , [0 Lim] , 'Color' , [.5 .5 .5] )
            set( gca , 'XLim' , [0  Lim] , 'YLim' , [0  Lim] )
           
            xlabel( 'Control' )
            ylabel( 'Observation' )
            grid on
        case 'T'
            Lim1 = min( min(Prc_obsr(1,:)) , min(Prc_ctrl(1,:)) ) ;
            Lim2 = max( max(Prc_obsr(end,:)) , max(Prc_ctrl(end,:)) ) ;
            H = plot( Prc_ctrl(:,1) , Prc_obsr(:,2) ) ;
%             opz_line = { '-' , ':' , '--' , '-.' } ;
%             for i = 1:I
%                 set( H(i) , 'LineWidth' , 1 , 'LineStyle' , opz_line{i} )
%             end
            hold on; plot( [Lim1 Lim2] , [Lim1 Lim2] , 'k' )
            set( gca , 'XLim' , [Lim1  Lim2] , 'YLim' , [Lim1  Lim2] )
            xlabel( 'Control' )
            ylabel( 'Observation' )
            grid on
    end
else
    H = [] ;
end

dwnsc_param.Prc_ctrl_bd = Prc_ctrl ;
dwnsc_param.Prc_obsr_bd = Prc_obsr ;
% % -----------------------------------------------------------------------
% % figure useful to understand the meaning of the following variables:
% Prc_C = sort( rand(99,1) -0.5 ) ;
% Prc_O = sort( rand(99,1) - 0.5) ;
% figure; plot( Prc_C , Prc_O , '.' )                                         % Q-Q plot
% hold on 
% plot( min(Prc_C) , min(Prc_O) , 'oc' )                                      % min
% plot( max(Prc_C) , max(Prc_O) , 'ok' )                                      % max
% plot( -1 , -1 + ( min(Prc_O) - min(Prc_C) ) , 'sc')                         % inf
% plot( 1 , 1 + ( max(Prc_O) - max(Prc_C) ) , 'sk')                           % sup
% plot( [-1 min(Prc_C)],[ -1 + ( min(Prc_O) - min(Prc_C) ) min(Prc_O)], '-c') % constant correction
% plot( [1 max(Prc_C)] , [1 + ( max(Prc_O) - max(Prc_C) ) max(Prc_O)] , '-k') % constant correction
% plot( [-1 1] , [-1 1] , 'r' )     
% legend( 'Q-Q plot' , 'min' , 'max' , 'inf' , 'sup' )
% xlabel('ctrl'); ylabel('obsr')
% % -----------------------------------------------------------------------

% Outside the table range (< 1 prc or > 99 prc), a constant correction is
% applied:
% e.g. if the 99th quantile of January temperature is corrected by +1°C,
% any temperature above this threshold is corrected by +1°C     
ctrl_min = min( Prc_ctrl )                       ;
ctrl_max = max( Prc_ctrl )                       ;
ctrl_inf = repmat ( C_I , 1 , size(Prc_ctrl,2) ) ;
ctrl_sup = repmat ( C_S , 1 , size(Prc_ctrl,2) ) ;
Prc_ctrl = [ ctrl_inf ; Prc_ctrl ; ctrl_sup ]    ;  % [ length(Prc)+2 , I ]

obsr_min = min( Prc_obsr )                    ;
obsr_max = max( Prc_obsr )                    ;
obsr_inf = ctrl_inf + ( obsr_min - ctrl_min ) ;     


switch var
    case 'P'
        obsr_inf( find( obsr_inf <0 ) ) = 0 ;       % avoid negative precipitation
end        
obsr_sup = ctrl_sup + ( obsr_max - ctrl_max ) ; 
Prc_obsr = [ obsr_inf ; Prc_obsr ; obsr_sup ] ;     % [ length(Prc)+2 , I ]


% A linear interpolation is applied between two percentiles: 
% the values of Prc_ctrl should be distinct, if they are not 
% (e.g. dry climate) an eps will be added.


ctrl_d = zeros(size(y_ctrl_original));

for i = 1:I
    
    
    Prc_ctrl( find( diff( Prc_ctrl(:,i) ) == 0  )+1 , i ) = Prc_ctrl( find( diff( Prc_ctrl(:,i) ) == 0  ) , i ) + eps .* find( diff( Prc_ctrl(:,i) ) == 0 ) ;
    % data to be modifies (x' = previous)        = data before the one to be modified                 + eps .* 99 (max, in case of 98 percentile  equal to 99 percentile, very difficult)
   
    [~, uidx] = unique(Prc_ctrl(:,i)); %  The grid vectors must contain unique points
    
    switch var
        
        case 'P'
            
            temp1=y_ctrl_original.*ID(:,i);
            temp1(temp1<TP)=0;
            temp1(temp1>=TP)=1;
            idx_temp1=logical(temp1);
            
        case 'T'
            
           idx_temp1=logical(ID(:,i)); 
    end
    
   ctrl_d(idx_temp1)= interp1( Prc_ctrl(uidx,i) , Prc_obsr(uidx,i) ,y_ctrl_original(idx_temp1)); 
  
   

end

% Vectors of downscaled time-series [ h , 1 ]:


dwnsc_param.Prc_ctrl = Prc_ctrl ;
dwnsc_param.Prc_obsr = Prc_obsr ;


end 