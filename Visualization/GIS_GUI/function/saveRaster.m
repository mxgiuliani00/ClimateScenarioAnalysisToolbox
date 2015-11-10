function  []  =  saveRaster (MatrData, Header, namefile)
% 
% description: this function is used to print the data in raster able to be
% read in ArcGIS, taking the processed matrix data and stored header as 
% inputs.


% Input
% MatrData  : processed matrix data
% Header    : stored info. of header from the original raster input file.
% namefile  : output file name


% Author      :  Yu Li
% Date        :  5th July, 2013
% Institution :  DEIB in Politecnico
% note: modified by MatteoG (15/7/2013)


[m, n] = size(MatrData);

fileID = fopen(namefile,'w');    

for i = 1: 6                % by assumption we have 6 lines of header
    fprintf(fileID,'%s %s \n',Header{1}{i}, Header{2}{i});
end                         % finish writing the header into file

for i = 1: m
    for j = 1: n
        if( isnan(MatrData(i,j)) )
            fprintf(fileID,'%d ',-9999);  % change the data type into the desidred one
        else
            %fprintf(fileID,'%5.2f ',MatrData(i,j));  % change the data type into the desidred one
            fprintf(fileID,'%g ',MatrData(i,j));  % change the data type into the desidred one
        end
    end
     fprintf(fileID,'\n');                  
end

fclose(fileID);


