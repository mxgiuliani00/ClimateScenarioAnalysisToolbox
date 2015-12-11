function  []  =  raster_print (MatrData, Header)
% 
% description: this function is used to print the data in raster able to be
% read in ArcGIS, taking the processed matrix data and stored header as 
% inputs.


% Input
% MatrData  : processed matrix data

% Header    : stored info. of header from the original raster input file.


% Author      :  Yu Li
% Date        :  5th July, 2013
% Institution :  DEIB in Politecnico


[m, n] = size(MatrData);

fileID = fopen('output.txt','a');    % change the name of file output in the future

for i = 1: 6                % by assumption we have 6 lines of header
    fprintf(fileID,'%s %s \n',Header{1}{i}, Header{2}{i});
end                         % finish writing the header into file

for i = 1: m
    for j = 1: n
     fprintf(fileID,'%5.2f',MatrData(i,j));  % change the data type into the desidred one
    end
     fprintf(fileID,'\n');                  
end

fclose(fileID);


