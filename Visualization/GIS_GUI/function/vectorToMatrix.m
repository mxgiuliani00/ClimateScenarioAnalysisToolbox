function M = vectorToMatrix( y, X )

% M = vectorToMatrix( y, X )
%
% This function converts the vector used in the model into a matrix with
% the same dimensions of the domain X, which can be then printed as a 
% raster file (the header with the geographical information must be added).
% Example for the Muzza model: the model works with vectors of 11667
% elements that are converted in a matrix of 182x202 cells.
%
% Input:    y = data vector
%           X = domain matrix (the one obtained from muzzagrid.asc)
%
% Output: 	M = data re-organized in a matrix according to the domain X
%
% MatteoG, 10/07/2013

l = length(y);
[r,c] = size(X) ;

if(l ~= sum( ~isnan(X(:)) ))
    error('Wrong input dimension: the length of y must be equal to the number of non-zero elements in x ');
end

% extraction of the elements in M according to the domain x
k = 1;
M = nan(r,c);
for i = 1:r
    for j = 1:c
        if X(i,j)==1
            M(i,j)  = y(k);
            k = k + 1;
        end
    end
end  



