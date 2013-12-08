function Y = cell2mat_pileup(X, mask)
% Y = cell2mat_pileup(X, mask)
% Convert a cell of matrices to a 3d array by piling up the 2d matrices and by applying the mask
% on each X{n}.

[M,N] = size(X{1});

L = length(X);

Y = zeros(M, N, L);
for n=1:L
    if nargin<2 || isempty(mask)
        Y(:,:,n) = X{n};
    else
        Y(:,:,n) = X{n}.*mask;
    end
end
