function Y = matrix_as_linop(A, X, tflag)
    if strcmp(tflag, 'notransp')
        Y = A*X;
    else
        Y = A'*X;
    end
end
