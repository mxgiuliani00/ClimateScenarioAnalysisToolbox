% convert vector to matrix
[filename_input, file_path]  = uigetfile('*.txt','Pick the file name to convert');
fielname_output = uiputfile('*.asc','Pick the file name to save as');

%% main execution part
file = [file_path filename_input];
data = load(file, '-ascii');


% read raster and split the header from the matrix
[muzza_domain, header]  =  loadRaster('muzzagrid.asc', 182); 
% muzza_domain is a matrix containing 1 in the cells belonging to the Muzza
% district and NaN outside.

M = vectorToMatrix( data, muzza_domain ) ;
% convert matrix in vector (11667 elements)

% save M as a raster
saveRaster(M, header, [fielname_output, '.asc']);

