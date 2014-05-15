function [Xc] = ckeep(X, N)
% Keep the central part (length of N) of an array X

M=length(X);

if N>=M
    Xc = X;
else
    idx = ceil(M/2) - floor(N/2);
    Xc = X(idx:(idx+N-1));
end