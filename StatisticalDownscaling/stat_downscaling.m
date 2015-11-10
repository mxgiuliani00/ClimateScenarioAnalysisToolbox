function [ scen_d, ctrl_d ] = stat_downscaling(obsr, ctrl, scen, varName, varargin)
%Statistical downscaling of climate change projections
%   [ scen_d, ctrl_d ] = stat_downscaling(obsr, ctrl, scen, varargin)
%   peroforms the statisical downscaling given the control time series and
%   the scenarios to be downscaled, and return both the downscaled control
%   and scenarios time-series. If otherwise specified, the default settings is to use
%   quantile-based method on the yearly temporal basis.
%
%   Input:
%
%     Name          Description
%
%     obsr          observed time-series (vector) with natural calendar.
%     ctrl          control time-series (vector).
%     scen          projected future scenario (vector).
%     varName       downscaling variables, 'P' for precipitation or 'T' for
%                   Temperature (C) ;
%     varargin{1}   optional parameter for the downscaling methods, default
%                   is quantile-based method ('qq').
%     varargin{2}   optional parameter for the temporal basis, default is
%                   on yearly basis ('annual').
%
%   Output:
%
%     Name          Description
%
%     scen_d        downscaled projected time-series (vector). If input
%                   'scen' is natural calendar, 'scen_d' will be 365-day
%                   calendar instead.
%     ctrl_d        downscaled control time-series (vector). If input
%                   'ctrl' is natural calendar, 'ctrl_d' will be 365-day
%                   calendar instead.
%
%   Examples:
%
%     Temperature-like time series (ctrl over estimate sistematically obsr)
%     noDays = 1096 ;
%     x = [1:length(noDays)]' ;
%     y_obs = sin(2*pi/365*x+143)*15+rand(size(x))   ;
%     y_bck = sin(2*pi/365*x+143)*15+rand(size(x))+5 ;
%
%     figure; plot( y_obs, '.-' )
%     hold on; plot( y_bck, '.-r' )
%     legend( 'obervation', 'backcast' )
%     xlabel('Days')
%     ylabel('Temperature [C]')
%
%     [ycdf_obs,xcdf_obs] = ecdf(y_obs) ;
%     figure; plot( xcdf_obs, ycdf_obs )
%     [ycdf_bck,xcdf_bck] = ecdf(y_bck) ;
%     hold on; plot( xcdf_bck, ycdf_bck ,'r' )
%     legend( 'obervation', 'backcast' )
%     xlabel('Temperature [C]')
%     ylabel('Empirical cdf')
%
%     [ y_bck_d, y_frc_d ] = stat_downscaling( y_obs , y_bck , y_frc , 'T' , 'qq', 'annual' ) ;
%
%     t_frc = [time_date2JD(01,01,2101):time_date2JD(31,12,2103)]' ;
%     t_dates_frc = time_JD2date(t_frc) ;
%     x_frc = [1:length(t_frc)]' ;
%     y_frc = sin(2*pi/365*x_frc+143)*15+2*rand(size(x_frc))+15 ;
%     [ y_frc_d ] = apply_downscaling_qq( y_frc , t_dates_frc , dwnsc_param.Prc_obsr , dwnsc_param.Prc_ctrl) ;
%
%     check: if you downscale y_bck --> qq plot must be a line 45 degrees
%     figure; plot( [y_obs y_bck_d] ); title('After downscaling')
%     figure; qqplot(y_bck_d,y_obs); grid on; title('After downscaling')

%   Author: Yu Li at Politecnico di Milano
%   Copyright 2015.

methodName = varargin{1} ;
tempBasis  = varargin{2} ;


switch tempBasis
  case 'annual'
    maskType = 1  ;
  case 'seasonal'
    maskType = 4  ;
  case 'monthly'
    maskType = 12 ;
  otherwise
    error(['Valid parameter name of temporal basis are:', ...
           '''annual'', ''seasonal'' or ''monthly'''])
end
%
% noYrs = floor(length(obsr)/365) ;
%
% if noYrs*365 ~= length(obsr)
%   obsr = obsr(1: noYrs*365) ;
%   obsr = reshape(obsr, 365, []) ;
% end
%
% tLength_ctrl = length(ctrl) ;
% tLength_scen = length(scen) ;
%
% if noYrs ~= floor(tLength_ctrl/365) && noYrs ~= floor(tLength_ctrl/360)
%   warning(['The No. of years between the observation and control time', ...
%     'series may be different'])
% end

% Trimming the input time-series to have consistent dimension
[ obsr, ctrl, scen ] = trim_timeseries( obsr, ctrl, scen )  ;

% Prepare the temporal mask
[nRows, ~] = size(ctrl) ;
if nRows == 360
  mask = 30*ones(12,1) ;
else
  mask = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]' ;
end

% Perform the downscaling

ctrStep = length(mask) / maskType ;

endIdx = cumsum(mask)               ;
begIdx = [ 1; endIdx(1:end-1) + 1 ] ;

ctrl_d = zeros(size(ctrl)) ;
scen_d = zeros(size(scen)) ;

for i = 1: maskType

  idx1 = begIdx( 1 + ctrStep*(i-1) ) ;
  idx2 = endIdx( ctrStep*i )         ;

  obsrTmp = obsr(idx1:idx2, :) ;
  ctrlTmp = ctrl(idx1:idx2, :) ;
  scenTmp = scen(idx1:idx2, :) ;

  switch methodName
    case 'qq'
      [ scen_d_tmp, ctrl_d_tmp ] = qq_downscaling( obsrTmp(:), ctrlTmp(:), scenTmp(:), varName ) ;
    otherwise
     error(['Valid downscaling methods are: ''qq'' '])  % TODO: update the methods
  end

  ctrl_d(idx1:idx2, :) = reshape(ctrl_d_tmp, idx2-idx1+1, []) ;
  scen_d(idx1:idx2, :) = reshape(scen_d_tmp, idx2-idx1+1, []) ;

end

ctrl_d = ctrl_d(:) ;
scen_d = scen_d(:) ;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subfunctions

%% Trimming the time-series of input variables
function [ obsr_new, ctrl_new, scen_new ] = trim_timeseries( obsr_old, ctrl_old, scen_old )

  noYrs_ctrl = floor(length(obsr_old)/365) ;

  tLength_obsr = length(obsr_old) ;
  tLength_ctrl = length(ctrl_old) ;
  tLength_scen = length(scen_old) ;

  if noYrs_ctrl*365 ~= tLength_obsr
    obsr_old = obsr_old(1: noYrs_ctrl*365) ;
    obsr_new = reshape(obsr_old, 365, [])  ; % [ 365 x No. of years ]
  end

  if noYrs_ctrl ~= floor(tLength_ctrl/365) && noYrs_ctrl ~= floor(tLength_ctrl/360)
    error(['The No. of years between the observation and control time', ...
             'series may be different'])
  end

  if mod(tLength_ctrl, 360) == 0                 % 360-day calendar
    ctrl_new = reshape(ctrl_old, 360, []) ;      % [ 360 x No. of years ]

    obsr_new = obsr_new(1:360, :) ;              % [ 360 x No. of years ]
  else                                           % 365-day calendar or natural calendar
    ctrl_old = ctrl_old(1:noYrs_ctrl*365) ;
    ctrl_new = reshape(ctrl_old, 365, []) ;      % [ 365 x No. of years ]
  end

  if mod(tLength_scen, 360) == 0                 % 360-day calendar
    scen_new = reshape(scen_old, 360, []) ;      % [ 360 x No. of years ]
  else                                           % 365-day calendar or natural calendar
    noYrs_scen = floor(length(scen_old)/365) ;
    scen_old   = scen_old(1:noYrs_scen*365)  ;
    scen_new   = reshape(scen_old, 365, [])  ;   % [ 365 x No. of years ]
  end

end


%% Quantile-based downscaling, adopted from DanielaA's script and has been
%  verified using the examples above

function [ scen_d, ctrl_d ] = qq_downscaling(obsr, ctrl, scen, varName)

  switch varName
    case 'P'
      C_I =   0 ; % mm/month
      C_S = 300 ; % mm/month
    case 'T'
      C_I = -50 ; % °C
      C_S =  50 ; % °C
  end

  Prc = [1:99]' ; % sampling percentile, to be improved with higher sampling frequency ?

  Prc_obsr = prctile(obsr, Prc) ;
  Prc_ctrl = prctile(ctrl, Prc) ;

  ctrl_min = min(Prc_ctrl)  ;
  ctrl_max = max(Prc_ctrl)  ;

  Prc_ctrl = [ C_I ; Prc_ctrl ; C_S ] ;  % [ length(Prc)+2 , 1 ]

  obsr_min = min(Prc_obsr)                 ;
  obsr_max = max(Prc_obsr)                 ;
  obsr_inf = C_I + ( obsr_min - ctrl_min ) ;

  if varName == 'P'
    obsr_inf(obsr_inf <0) = 0 ;       % avoid negative precipitation
  end

  obsr_sup = C_S + ( obsr_max - ctrl_max ) ;
  Prc_obsr = [ obsr_inf ; Prc_obsr ; obsr_sup ] ;     % [ length(Prc)+2 , 1 ]

  % A linear interpolation is applied between two percentiles:
  % the values of Prc_ctrl should be distinct.

  [ Prc_ctrl, idxTmp ] = unique(Prc_ctrl) ;
  Prc_obsr           = Prc_obsr(idxTmp) ;

  ctrl_d = interp1( Prc_ctrl , Prc_obsr , ctrl ) ;
  scen_d = interp1( Prc_ctrl , Prc_obsr , scen ) ;

end

%% TODO: add other methods
