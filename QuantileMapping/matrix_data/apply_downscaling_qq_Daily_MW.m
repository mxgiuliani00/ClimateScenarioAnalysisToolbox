function [ scen_d ] = apply_downscaling_qq_Daily_MW( scen ,var, Prc_obsr , Prc_ctrl, TP)

% [ scen_d ] = apply_downscaling_qq( scen, tdates_scen, Prc_obsr, Prc_ctrl)
%
% The function downscales precipitation or temperature timeseries simulated
% through a climate model. Downscaling is made by the Quantile Method. 
% 
% The downscaling function is given by a Q-Q plot estimated on a daily (365 
% correction functions) computed on a moving window. The window has an
% amplitude 2 x f and is centered in each day of the year.
% A Q-Q plot is a probability plot, a kind of graphical method for 
% comparing two probability distributions, by plotting their quantiles
% against each other. 
% The correction function is based on the 99 percentiles of the 
% distribution of the observation and of the control timeserie.
% A linear interpolation is applied between two percentiles.
% Outside that range (i.e. < 1 prc or > 99 prc), a constant 
% correction is applied: if the 99th quantile of temperature is corrected
% by +1°C, any temperature above this threshold is corrected by +1°C. 
% 
%
% Input:
%   - scen        = scenario time-serie                          - vector [ h_scen , 1 ] 
%   - var         = type of variable to be downscaled.           - string
%                   It can assume two values only: 'P' = precipitation
%                                                  'T' = temperature
%   - Prc_obsr    = percentiles of the observation               - matrix [101,I]
%   - Prc_ctrl    = percentiles of the control timeserie         - matrix [101,I]
%                   
% Output:
%   - scen_d      = downscaled scenario timeserie                - vector [ h_scen , 1 ]
% 




% matrix  [ h , 365 ] of logical value: 1s in the j-th column stand for
% data of the time serie belonging to the j-th day of the year. The matrix
% ID is the index of the day of the year for which the correction function
% is estimated. 

I = size(Prc_ctrl,2) ;
Ny = floor( length(scen)/365 ) ;
ID = repmat( eye(365), Ny, 1 ) ;
scen_d = zeros(size(scen));

% A linear interpolation is applied between two percentiles: 
% the values of Prc_ctrl should be distinct, if they are not 
% (e.g. dry climate) an eps will be added.
for i = 1:I
     [~, uidx] = unique(Prc_ctrl(:,i)); 
      
     switch var
        
        case 'P'
            
 % In case of P, data for dirzzling need to be corrected. Precipitation
 % below 1 mm/day is set to zero and the rest above the threshold is used
 % to estime percentile.
            
            temp1=scen.*ID(:,i);
            temp1(temp1<TP)=0;
            temp1(temp1>=TP)=1;
            idx_temp1=logical(temp1);
            
        case 'T'
            
           idx_temp1=logical(ID(:,i)); 
           
     end
%     scen_d( find(ID(:,i)) ) = interp1( Prc_ctrl(uidx,i) , Prc_obsr(uidx,i) , scen( find(ID(:,i)) ) ) ;
      scen_d( idx_temp1 ) = interp1( Prc_ctrl(uidx,i) , Prc_obsr(uidx,i) , scen( idx_temp1 ) ) ;
end

% Vectors of downscaled time-series [ h , 1 ]:
scen_d = scen_d' ;
