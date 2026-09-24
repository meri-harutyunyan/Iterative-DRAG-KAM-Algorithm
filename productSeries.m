function seriesAB = productSeries(seriesA, seriesB)
%-------------------------------------------------------------------------
% ProductSeries: Multiplies two operator series up to a given order.
%
%
%   Given the series of operator A and B, produces the series of their product. 
%   A series is defined as:
%       X = Σ_{i=1}^∞ ε^(i-1) * seriesX(:,:,i)
%   The product serie will be:
%       seriesAB(:,:,i) = Σ_{j=1}^i seriesA(:,:,j) * seriesB(:,:,i+1-j)
% 
% 
%   INPUTS:
%       seriesA   - First operator series,  size [3, 3, Na].
%       seriesB   - Second operator series, size [3, 3, Nb].
%
%   OUTPUT:
%       seriesAB  - Product operator series, size [3, 3, maxOrder].
%------------------------------------------------------------------------- 


% Ensure maxOrder does not exceed input series lengths
maxOrder = min([size(seriesA,3), size(seriesB,3)]);
numLevels = size(seriesA,2);
% Compute seriesAB
seriesAB = zeros(numLevels, numLevels, maxOrder);
for i = 1:maxOrder
    for j = 1:i
        seriesAB(:,:,i) = seriesAB(:,:,i) + seriesA(:,:,j) * seriesB(:,:,i+1-j);
    end
end

end
