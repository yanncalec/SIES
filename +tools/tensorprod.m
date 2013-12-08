function XY = tensorprod(X, Y)
% XY = tensorprod(X, Y)
% Calculate the tensor product between two vectors X and Y. The output is a
% matrix with (m,n)-th coefficient given by X(m)*Y(n)


    X = reshape(X, [], 1); Y = reshape(Y, 1, []);
    XY = repmat(X, 1, length(Y)) .* repmat(Y, length(X), 1);
    
    % % Or by mex function:
    % XY = tools.tensorprod_mex(X, Y);

end

