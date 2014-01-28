function [A, D] = wavedec2(X, Nx1, Ny1, hcoeff, hrange, nlevel)
% 2D wavelet decomposition with explicit position index range and mirror filter index range. Unlike
% the discrete wavelet transform in Matlab or WaveLab, the position indexes are tracked during the
% transform. This is particularly useful when the input X is the value of a function evaluated on a
% 2d grid and when we need to know explictly the position index of the output wavelet coefficients.
%
% INPUTS:
% X: 2D matrix to be decomposed
% Nx1, Ny1: starting index 
% hcoeff: qmf coefficients
% hrange: index range of the qmf
% nlevel: level of decomposition
% OUTPUTS:
% A,D: approximation and detail coefficients. D{end,1}..D{end,3} are the finest
% scale (Jmin), while D{1,1}..D{1,3} are the coarsest scale.
    
    if nlevel==0
        W = X(:);
        A.coeff = X;
        A.rangex = Nx1;
        A.rangey = Ny1;
        D={};
        return;
    end
        
    [nrow, ncol] = size(X);

    % No wavelet transform if the length of signal is smaller than that of qmf
    max_nlevel = floor(log2(min(nrow, ncol)/length(hcoeff)));

    if nrow==1 || ncol==1
        error('Input matrix must have 2 dimension.');
    elseif nlevel > max_nlevel
        error('Decomposition level is too large.');
    end

    % High frequency reconstruction filter
    gcoeff = wavelet.OrthoWvl.qmf_h2g(hcoeff, hrange);

    Acoeff = X; 
    D = cell(nlevel, 3);

    for j=nlevel:-1:1
        [A, D{j,1}, D{j,2}, D{j,3}] = wavelet.OrthoWvl.dwt2(Acoeff, Nx1, Ny1, hcoeff, hrange, ...
                                                          gcoeff); 
        Nx1 = A.rangex(1);
        Ny1 = A.rangey(1);
        
        Acoeff = A.coeff;
    end
end

% function parity = check_parity(A, B)
%     if mod(A,2) && mod(B,2)
%         parity = 'oo';
%     elseif mod(A,2) && ~mod(B,2)    
%         parity = 'oe';
%     elseif ~mod(A,2) && mod(B,2)
%         parity = 'eo';
%     else
%         parity = 'ee';
%     end
% end

