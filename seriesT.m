function T = seriesT(M,WM,numOrders)
% -------------------------------------------------------------------------
%seriesT: Construct operator series T from powers of WM.
%
%   Builds the operator series expansion of T using powers of WM:
%       T = exp(ε^M WM)
%       T(:,:,i) = (1 / j!) * WM^j
%   where j = (i-1)/M and i runs in steps of M.
%
%   INPUTS:
%       M         - Iteration number.
%       WM        - W at M-th iteration, 3x3 matrix.
%       numOrders - Maximum expansion order.
%
%   OUTPUT:
%       T         - Transformation operator series, size 3x3xnumOrders.
%
% -------------------------------------------------------------------------

numLevels = size(WM,2);
T = zeros(numLevels,numLevels,numOrders);
for i=1:M:numOrders
    j = (i-1)/M;
    T(:,:,i) = WM^(j)/factorial(j);
end

end