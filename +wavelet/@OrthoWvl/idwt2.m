function out = idwt2(A, D1, D2, D3, hcoeff, hrange, gcoeff)
    Nh1 = hrange(1); Nh2 = hrange(2);
    LoR = hcoeff;
    HiR = gcoeff;

    %% Apprx coefficients (Lo in row, Lo in column)
    toto = dyadup(A.coeff);

    Acoeff = conv2(LoR, LoR, toto, 'full'); 
    Nx1 = A.rangex(1); Ny1 = A.rangey(1);

    [Px_A, Py_A] = wavelet.OrthoWvl.matrix_idx_range(Acoeff, Nh1+2*Nx1, Nh1+2*Ny1);
    
    %% D1 coefficients (Hi in row, Lo in column)
    toto = dyadup(D1.coeff);

    D1coeff = conv2(HiR, LoR, toto, 'full');
    Nx1 = D1.rangex(1); Ny1 = D1.rangey(1);

    [Px_D1, Py_D1] = wavelet.OrthoWvl.matrix_idx_range(D1coeff, Nh1+2*Nx1, 1-Nh2+2*Ny1);

    %% D2 coefficients (Lo in row, Hi in column)
    toto = dyadup(D2.coeff);

    D2coeff = conv2(LoR, HiR, toto, 'full');
    Nx1 = D2.rangex(1); Ny1 = D2.rangey(1);

    [Px_D2, Py_D2] = wavelet.OrthoWvl.matrix_idx_range(D2coeff, 1-Nh2+2*Nx1, Nh1+2*Ny1);

    %% D3 coefficients (Hi in row, Hi in column)
    toto = dyadup(D3.coeff);

    D3coeff = conv2(HiR, HiR, toto, 'full');
    Nx1 = D3.rangex(1); Ny1 = D3.rangey(1);

    [Px_D3, Py_D3] = wavelet.OrthoWvl.matrix_idx_range(D3coeff, 1-Nh2+2*Nx1, 1-Nh2+2*Ny1);
    
    %% Sum of A, D1, D2, D3
    col0 = min([Px_A(1), Px_D1(1), Px_D2(1), Px_D3(1)]);
    col1 = max([Px_A(2), Px_D1(2), Px_D2(2), Px_D3(2)]);

    row0 = min([Py_A(1), Py_D1(1), Py_D2(1), Py_D3(1)]);
    row1 = max([Py_A(2), Py_D1(2), Py_D2(2), Py_D3(2)]);
    
    ndim = [row1-row0+1, col1-col0+1];
    Ae = tools.extend_matrix(Acoeff, ndim, Py_A(1)-row0+1, Px_A(1)-col0+1); clear Acoeff;
    D1e = tools.extend_matrix(D1coeff, ndim, Py_D1(1)-row0+1, Px_D1(1)-col0+1);  clear D1coeff;
    D2e = tools.extend_matrix(D2coeff, ndim, Py_D2(1)-row0+1, Px_D2(1)-col0+1);  clear D2coeff;
    D3e = tools.extend_matrix(D3coeff, ndim, Py_D3(1)-row0+1, Px_D3(1)-col0+1);  clear D3coeff;

    % clear Acoeff D1coeff D2coeff D3coeff;
    
    out.coeff = Ae+D1e+D2e+D3e;
    out.rangex = [col0, col1];
    out.rangey = [row0, row1];
end

function toto = dyadup(X)
    [nrow, ncol] = size(X);

    toto = zeros(2*nrow-1, 2*ncol-1);
    toto(1:2:end, 1:2:end) = X;

    % [Px,Py] = wavelet.OrthoWvl.matrix_idx_range(toto, 2*Nx, 2*Ny);

    % out.coeff = coeff; % down-sampled coefficient
    % out.rangex = Px; % position index in column direction of coeff
    % out.rangey = Py; % position index in row direction of coeff
end
