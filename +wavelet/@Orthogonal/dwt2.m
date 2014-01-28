function [A, D1, D2, D3] = dwt2(X, Nx1, Ny1, hcoeff, hrange, gcoeff)    
    [nrow, ncol] = size(X);

    Nh1 = hrange(1); Nh2 = hrange(2);

    LoD = hcoeff(end:-1:1);    
    HiD = gcoeff(end:-1:1);

    % conv2(hrow, hcol, X) convolves first vertically (in row direction) X then
    % convolves horizontally (in column direction).

    % Some rules:
    % If it is low-pass filtering in x(or y), then Ncol=Nx1-Nh2, otherwise Ncol=Nx1+Nh1,
    % Same thing for Nrow.
    
    %% Apprx coefficients (Lo in row, Lo in column)
    toto = conv2(LoD, LoD, X, 'full'); 
    
    % Dyadic downsampling
    Nrow = Ny1-Nh2; Ncol = Nx1-Nh2; % <--Danger 1

    pstr = [num2str(mod(Nrow, 2)), num2str(mod(Ncol, 2))]; % row-column
    switch pstr
      case '11' % 'row-odd, col-odd'
        coeff = toto(2:2:end, 2:2:end); % <--Danger 2
        [Px, Py] = wavelet.OrthoWvl.matrix_idx_range(coeff, (Ncol+1)/2, (Nrow+1)/2); % <--Danger 3
      case '10' % 'row-odd, col-even'
        coeff = toto(2:2:end, 1:2:end);
        [Px, Py] = wavelet.OrthoWvl.matrix_idx_range(coeff, (Ncol)/2, (Nrow+1)/2);
      case '01' % 'row-even, col-odd'
        coeff = toto(1:2:end, 2:2:end);
        [Px, Py] = wavelet.OrthoWvl.matrix_idx_range(coeff, (Ncol+1)/2, (Nrow)/2);        
      case '00' % 'row-even, col-even'
        coeff = toto(1:2:end, 1:2:end);
        [Px, Py] = wavelet.OrthoWvl.matrix_idx_range(coeff, (Ncol)/2, (Nrow)/2);
    end
    
    A.coeff = coeff; A.rangex = Px; A.rangey = Py;
    
    %% D1 coefficients (Hi in row, Lo in column)
    toto = conv2(HiD, LoD, X, 'full'); 
    
    % Dyadic downsampling
    Nrow = Ny1+Nh1; Ncol = Nx1-Nh2; %<--Danger 1

    pstr = [num2str(mod(Nrow, 2)), num2str(mod(Ncol, 2))]; % row-column
    switch pstr
      case '11' % 'row-odd, col-odd'
        coeff = toto(1:2:end, 2:2:end); %<--Danger 2
        [Px, Py] = wavelet.OrthoWvl.matrix_idx_range(coeff, (Ncol+1)/2, (Nrow-1)/2); %<--Danger 3
      case '10' % 'row-odd, col-even'
        coeff = toto(1:2:end, 1:2:end);
        [Px, Py] = wavelet.OrthoWvl.matrix_idx_range(coeff, (Ncol)/2, (Nrow-1)/2);
      case '01' % 'row-even, col-odd'
        coeff = toto(2:2:end, 2:2:end);
        [Px, Py] = wavelet.OrthoWvl.matrix_idx_range(coeff, (Ncol+1)/2, (Nrow)/2);        
      case '00' % 'row-even, col-even'
        coeff = toto(2:2:end, 1:2:end);
        [Px, Py] = wavelet.OrthoWvl.matrix_idx_range(coeff, (Ncol)/2, (Nrow)/2);
    end

    D1.coeff = coeff; D1.rangex = Px; D1.rangey = Py;

    %% D2 coefficients (Lo in row, Hi in column)
    toto = conv2(LoD, HiD, X, 'full'); 
    
    % Dyadic downsampling
    Nrow = Ny1-Nh2; Ncol = Nx1+Nh1; %<--Danger 1

    pstr = [num2str(mod(Nrow, 2)), num2str(mod(Ncol, 2))]; % row-column
    switch pstr
      case '11' % 'row-odd, col-odd'
        coeff = toto(2:2:end, 1:2:end); %<--Danger 2
        [Px, Py] = wavelet.OrthoWvl.matrix_idx_range(coeff, (Ncol-1)/2, (Nrow+1)/2); %<--Danger 3
      case '10' % 'row-odd, col-even'
        coeff = toto(2:2:end, 2:2:end);
        [Px, Py] = wavelet.OrthoWvl.matrix_idx_range(coeff, (Ncol)/2, (Nrow+1)/2);
      case '01' % 'row-even, col-odd'
        coeff = toto(1:2:end, 1:2:end);
        [Px, Py] = wavelet.OrthoWvl.matrix_idx_range(coeff, (Ncol-1)/2, (Nrow)/2);        
      case '00' % 'row-even, col-even'
        coeff = toto(1:2:end, 2:2:end);
        [Px, Py] = wavelet.OrthoWvl.matrix_idx_range(coeff, (Ncol)/2, (Nrow)/2);
    end

    D2.coeff = coeff; D2.rangex = Px; D2.rangey = Py;

    %% D3 coefficients (Hi in row, Hi in column)
    toto = conv2(HiD, HiD, X, 'full'); 
    
    % Dyadic downsampling
    Nrow = Ny1+Nh1; Ncol = Nx1+Nh1; %<--Danger 1

    pstr = [num2str(mod(Nrow, 2)), num2str(mod(Ncol, 2))]; % row-column
    switch pstr
      case '11' % 'row-odd, col-odd'
        coeff = toto(1:2:end, 1:2:end); %<--Danger 2
        [Px, Py] = wavelet.OrthoWvl.matrix_idx_range(coeff, (Ncol-1)/2, (Nrow-1)/2); %<--Danger 3
      case '10' % 'row-odd, col-even'
        coeff = toto(1:2:end, 2:2:end);
        [Px, Py] = wavelet.OrthoWvl.matrix_idx_range(coeff, (Ncol)/2, (Nrow-1)/2);
      case '01' % 'row-even, col-odd'
        coeff = toto(2:2:end, 1:2:end);
        [Px, Py] = wavelet.OrthoWvl.matrix_idx_range(coeff, (Ncol-1)/2, (Nrow)/2);        
      case '00' % 'row-even, col-even'
        coeff = toto(2:2:end, 2:2:end);
        [Px, Py] = wavelet.OrthoWvl.matrix_idx_range(coeff, (Ncol)/2, (Nrow)/2);
    end
    D3.coeff = coeff; D3.rangex = Px; D3.rangey = Py;
end

