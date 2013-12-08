function [Xa, err, idx] =  Best_Napprx(X, N)
% [Xa, err, idx] =  Best_Napprx(X, N)
% Best N term approximation of a vector
% Inputs:
% X, N: the vector and the number of approximation terms
% Outputs:
% Xa: vector of approximation
% err: error of approximation
% idx: index of approximation

[Xs, idx] = sort(abs(X(:)), 'descend');
Xa = X;

if N<numel(Xa) && N>0
    Xa(idx(N+1:end)) = 0;
end

err = norm(Xa - X, 'fro');
idx = idx(1:N);