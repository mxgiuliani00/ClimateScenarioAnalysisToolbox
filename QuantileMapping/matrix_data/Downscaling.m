
% This script shows how to perform a time-varying quantile-quantile
% downscaling by correcting the outputs of a climate model to match the
% statistical characteristics of the local observations. Beside being
% quantile-specific, the correction function varies across days and is
% computed using a seasonal moving window (the dimension of the moving
% window can be adjusted).
%
% INPUTS
% - observation: vector/matrix of observations(if different stations rows:days x columnS:Nstation) same length of control
% - control: vector/matrix of RCM output for control(if different stations rows:days x columnS:Nstation) same length of obs
% - future: vector/matrix of RCM output for future(if different stations rows:days x columnS:Nstation)
% - var: 'P' for precipitation, 'T' for temperature
%
% OUTPUTS
% - ControlD: vector/matrix of downscaled climate model output over
%               the control period
% - FutureD: vector/matrix of downscaled future scenarios
%
% FUNCTIONS NEEDED:
% - compute_downscaling_qq_Daily_MW
% - apply_downscaling_qq_Daily_MW
%
% USER SETTINGS
% - f : semi-amplitude of the moving window (default is 45 days -> seasonal window)
% - TP: precipitation threshold (default is 1 mm/day)
%
% Copyright 2020 Environmental Intelligence Lab - Politecnico di Milano
% 
% Developers: Matteo Giuliani, Patrizia Zamberletti, Daniela Anghileri.

% This file is part of ClimateScenarioAnalysisToolbox repository.
% 
%     MATLAB_IterativeInputSelection_with_RTree-c is free software: you can redistribute 
%     it and/or modify it under the terms of the GNU General Public License 
%     as published by the Free Software Foundation, either version 3 of the 
%     License, or (at your option) any later version.     
% 
%     This code is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with MATLAB_IterativeInputSelection_with_RTree-c.  
%     If not, see <http://www.gnu.org/licenses/>.

%% Set workspace
clear
clc

%% Prepare input data and downscaling settings

% load data (observation, control period, future scenario)
load -ascii Como_Tdataset.txt   % 5-years trajectory of daily temperature
observation = Como_Tdataset(:,1);
control = Como_Tdataset(:,2);
future = Como_Tdataset(:,3);
var = 'T';

% define downscaling parameters
w = 45;     % semi-aplitude of the moving window [days]
TP = 1;     % precipitation treshold [mm/d]
 
%% Run dowscanling

VariableSpec = "%s";
switch var 
    case 'P'
        A1 = ' Precipitation [mm/d] ';
    case 'T'
        A1 = ' Temperature [C°] ';     
end

T= 365;

for ss=1:size(observation,2) %--> ss= #stations
    
    % data for station ss
    y_obs = observation(:,ss) ;
    y_ctrl = control(:,ss) ;
    
    % plot data
    figure; plot( y_obs, '.-k' )
    hold on; plot( y_ctrl, '.-r' )
    legend( 'observation', 'control' )
    xlabel('Time [days]')
    ylabel( sprintf(VariableSpec,A1))
    
    % CFD
    [ycdf_obs,xcdf_obs] = ecdf(y_obs) ;
    figure; plot( xcdf_obs, ycdf_obs, 'k' )
    [ycdf_ctrl,xcdf_ctrl] = ecdf(y_ctrl) ;
    hold on; plot( xcdf_ctrl, ycdf_ctrl ,'r' )
    legend( 'observation', 'control' )
    xlabel( sprintf(VariableSpec,A1) )
    ylabel('Empirical cdf')
    
   % DOWNSCALING - calibration phase
    [ y_ctrl_d, dwnsc_param , H ] = compute_downscaling_qq_Daily_MW(y_obs , y_ctrl , var , w , 0 ,TP) ;
    ControlD(:,ss)=y_ctrl_d;
    
    % check1: ctrl trajectory after downscaling
    figure; plot( y_obs, '.-k' )
    hold on; plot( y_ctrl, '.-r' )
    hold on; plot( y_ctrl_d, 'x--c' )
    legend( 'observation', 'control', 'control after downscaling' )
    xlabel('Days')
    ylabel(A1)
    
    % check2: if you downscale y_bck --> qq plot must be a line 45 degrees
    pc = prctile( y_ctrl_d(y_ctrl_d>TP) , [1:99] ) ;
    po = prctile( y_obs(y_obs>TP) , [1:99] ) ;
    Lim = max( max(po) , max(pc) ) ;
    figure; plot( pc,po,'x')
    hold on; plot( [0 Lim] , [0 Lim] , 'Color' , [.5 .5 .5] )
    
    % DOWNSCALING - projection phase
    y_scen = future(:,ss);
    [ y_scen_d ] = apply_downscaling_qq_Daily_MW( y_scen ,var, dwnsc_param.Prc_obsr , dwnsc_param.Prc_ctrl, TP) ;
    
    % check: scen trajectory after downscaling
    figure; plot( y_scen, '.-b' )
    hold on; plot( y_scen_d, 'x--m' )
    legend( 'scenario', 'scenario after downscaling' )
    xlabel('Days')
    ylabel(A1)
    
    FutureD(:,ss) = y_scen_d' ;    
end


