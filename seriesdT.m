function dT=seriesdT(M,WM,dWM, numOrders)
% -------------------------------------------------------------------------
% seriesdT: Construct operator series T from powers of WM.
%
%   Builds the operator series expansion of dT using powers of WM:
%       dT(:,:,i) = (1 / (j+1)!) * sum_{k=0}^{j} WM^k * dWM * WM^(j-k),
%   because generally WM and dWM don't cummute with each other.
%
%   INPUTS:
%       M         - Iteration number.
%       WM        - W at M-th iteration, 3x3 matrix.
%       dWM       - Time derivative of WM, 3x3 matrix.
%       numOrders - Maximum expansion order.
%
%   OUTPUT:
%       dT        - Transformation operator derivative series, size 3x3xnumOrders.
% -------------------------------------------------------------------------

numLevels = size(WM,2);
dT = zeros(numLevels,numLevels,numOrders);
for i=M+1:M:numOrders
    j = (i-M-1)/M;
    for k=0:j
         dT(:,:,i) = dT(:,:,i) + (WM^(k))*dWM*(WM^(j-k)) ;
    end
    dT(:,:,i) = dT(:,:,i)/factorial(j+1);
end

end