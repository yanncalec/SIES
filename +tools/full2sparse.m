function S = full2sparse(X, epsilon)
% S = full2sparse(X, epsilon)
% Convert a full array X to sparse form by truncating the coefficients of X
% smaller than (absolute value) max(X) * epsilon

    idx = abs(X) > (max(abs(X(:))) * epsilon);
    S = sparse(X .* idx);

end

