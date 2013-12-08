function S = l0norm(X)
    S = sum(abs(X(:))>0);
end