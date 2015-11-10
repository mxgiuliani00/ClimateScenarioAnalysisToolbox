function y = matrixToVector( M, X )

% y = matrixToVector( M, X )
%
% This function converts the matrix extracted from a raster (without the
% header) in the vector used in the model according to the domain X. 
% Example for the Muzza model: the raster file has 182x202 cells and is
% then converted in a column vector of 11667 elements.
%
% Input:    M = generic data matrix 
%           X = domain matrix (the one obtained from muzzagrid.asc)
%
% Output: 	y = data re-organized in a vector according to the
% implementation of the model
%
% MatteoG, 10/07/2013

[r,c] = size(M) ;
[r1,c1] = size(X) ;

if(r ~= r1 || c ~= c1)
    error('Wrong input dimension: M and x must have the same size');
end

% extraction of the elements in M according to the domain x
k = 1;
x = ~isnan(X);
y = nan(sum(x(:)),1);
for i = 1:r
    for j = 1:c
        if x(i,j)
            y(k)  = M(i,j);
            k = k + 1;
        end
    end
end  

end

