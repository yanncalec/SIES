function Y = extend_matrix(X, dim, row1, col1)
% Y = extend_matrix(X, dim, row1, col1)
% Extend a matrix X to a larer dimension 
% INPUTS:
% X: matrix to be extended
% dim: dimension to be extened to
% row1,col1: starting index (row and column) for X in the extended matrix
% If row1, col1 are not given, then the extension is symmetric by placing X in middle.
    
    [nrow, ncol] = size(X);

    if nargin < 3
        if dim(1)< nrow || dim(2) < ncol
            error('Dimension error');
        end

        Y = zeros(dim);
        r0 = floor(dim(1)/2)+1;
        c0 = floor(dim(2)/2)+1;

        r1 = r0 - floor(nrow/2);
        c1 = c0 - floor(ncol/2);

        Y(r1:r1+nrow-1, c1:c1+ncol-1) = X;        
    else        
        if row1+nrow-1>dim(1) || col1+ncol-1>dim(2)
            error('Dimension error');
        end

        Y = zeros(dim);
        Y(row1:row1+nrow-1, col1:col1+ncol-1) = X;
    end
end