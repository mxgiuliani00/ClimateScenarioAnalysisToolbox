function [rasterMatr, header]  =  raster_process(rasterFile, M)
% description: this function is used to process the raster output from
% ArcGIS software, and generates the numeric values in matrix fomrm and the
% header part.

% Input 
% rasterFile: the raster data from ArcGIS software
% M         : the No. of columns of the data contained in raster file.

% Output
% rasterMatr: matrix that containing the numeric values from raster file.
%             Nan values are set to zero;
% header    : array that contains the string information of header from
%             raster file.

% Author      :  Yu Li
% Date        :  5th July, 2013
% Institution :  DEIB in Politecnico

rasterMatr = [];
header     = {[]};

fileID     = fopen(rasterFile);
formatSpec = '%s %s';
N          = 6;        % by default, ArcGIS raster output has header of 6 lines  
% M = 462 ;            % 462 columns in the file

% subtract the header
header = textscan(fileID, formatSpec, N);

% subtract the numerical data;
while ~feof(fileID)
   rasterData          = textscan(fileID, '%f', M);   % in cell
   if (isempty(rasterData{1}))
       break
   end
   temp                = rasterData{1};
   rasterMatr(end+1,:) = temp';              
   %    process_data(rasterData);  % process_data is a function which
   %                                 process the 'rasterData ' into desired values
end
fclose(fileID);

rasterMatr(rasterMatr == -9999) = 0; % supress the NaN value into zeros


