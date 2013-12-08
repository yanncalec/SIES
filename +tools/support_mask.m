function M = support_mask(X)
% M = support_mask(X)
% Construct the 0-1 mask for the support of X
    
    idx = find(abs(X)>0);
    M = zeros(size(X));
    M(idx) = 1;
end