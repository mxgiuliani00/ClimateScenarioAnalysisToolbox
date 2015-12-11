function xx = cutNetCDFdata( x, mLon, mLat )

% xx = cutNetCDFdata( x, mLon, mLat )
%
% This function cuts the matrix stored in the matlab structure x (created
% via extractNetCDFdata function) for a given spatial domain.
%
% Input:    - x = matlab structure obtained from a NetCDF file using the
%                   extractNetCDFdata function 
%           - mLon = longitude domain (vector [minLon, maxLon])
%           - mLat = latitude domain (vector [minLat, maxLat])
% Output:   - xx = matlab structure containing the coordinates of the
%                   spatial domain of itnerests in .lon and .lat, along
%                   with the extracted time series in .value
%
% Last Update: MatteoG, 11/12/2015

% extracting domain from ncdf file
lon = x.lon;
lat = x.lat;

% creating domain mask
idxLon = ( lon > mLon(1) )&( lon < mLon(2) );
idxLat = ( lat > mLat(1) )&( lat < mLat(2) );
idxDomain = idxLon & idxLat;
idxDomain3D = repmat(idxDomain,[1,1,size(x.value,3)]);

% cutting netcdf
xx = x;
xx.lon = xx.lon(idxDomain) ;
xx.lat = xx.lat(idxDomain) ;
xx.value = xx.value(idxDomain3D) ;

end

% Copyright 2015 Yu Li and Matteo Giuliani, Politecnico di Milano
% M. Giuliani: matteo.giuliani@polimi.it - http://giuliani.faculty.polimi.it