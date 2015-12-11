function xx = extractNetCDFdata( filename, var )

% xx = extractNetCDFdata( filename, var )
%
% This function constructs a matlab structure containing the information
% from a NetCDF file.
%
% Input:    - filename = name of the NetCDF file to load (string)
%           - var = type of variable for conversion of unit of measures (string)
%                   the function currently works only for precipitation (pr) 
%                   and temperature (tas or ts)
% Output:   - xx = 5-fields matlab structure containing the NetCDF metadata
%                   in .headerInfo, longitude and latitude matrixes in .lon
%                   and .lat, temperature/precipitation matrix in .value,
%                   type of experiment (historical vs RCP) in .experiment
%
% Last Update: MatteoG, 11/12/2015

% collect the header information
fInfo = ncinfo(filename) ;
headerInfo = cell(length(fInfo.Attributes), 2) ;
for i = 1: length(fInfo.Attributes)
    headerInfo{i, 1} = fInfo.Attributes(1,i).Name  ;
    headerInfo{i, 2} = fInfo.Attributes(1,i).Value ;
end
xx.headerInfo = headerInfo ;

% read the coordinates
varLon = ncread(filename, 'lon') ;
varLat = ncread(filename, 'lat') ;

xx.lon = varLon ;
xx.lat = varLat ;

% read variables (surface temperature)
if strcmp(var,'tas') || strcmp(var,'ts')
    tas = ncread(filename, var) ; % °K
    xx.value = tas- 273.15 ; % °C
elseif strcmp(var,'pr')
    pr = ncread(filename, var) ; % kg/(m2 s) = mm/s
    xx.value = pr*3600*24 ; % mm/day
else 
    error('the selected variable is neither precipitation (pr) or temperature (tas,ts)')
end

% add metadata
xx.experiment = xx.headerInfo{5,2};

end

% Copyright 2015 Yu Li and Matteo Giuliani, Politecnico di Milano
% M. Giuliani: matteo.giuliani@polimi.it - http://giuliani.faculty.polimi.it