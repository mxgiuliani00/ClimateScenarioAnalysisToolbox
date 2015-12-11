
% script to test input-output functions that connect the model to raster
% files.
%
% MatteoG 15/07/2013

clear
clc

% read raster and split the header from the matrix
[muzza_domain, header]  =  loadRaster('muzzagrid.asc', 182); 
% muzza_domain is a matrix containing 1 in the cells belonging to the Muzza
% district and NaN outside.

% convert matrix in vector (11667 elements)
y = matrixToVector( muzza_domain, muzza_domain );

% convert vector to matrix
load -ascii comizi.txt
M = vectorToMatrix( comizi, muzza_domain ) ;

% save M as a raster
saveRaster(M, header, 'testFile.asc');

% clear
% clc
% 
% file = dir('crop*.txt');
% [muzza_domain, header]  =  loadRaster('muzzagrid.asc', 182);
% 
% for i = 1: length(file)
%     crop_choice = zeros(11667,1);
%     crop_choice = load(file(i).name);
%     M           = vectorToMatrix(crop_choice, muzza_domain ) ;
%     filename    = strcat(strtok(file(i).name,'.'),{'.asc'});
%     saveRaster(M, header, filename{1});
% end
