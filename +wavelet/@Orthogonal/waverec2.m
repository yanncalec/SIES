function A0 = waverec2(A, D, hcoeff, hrange)
% High frequency decomposition filter
    
    gcoeff = wavelet.OrthoWvl.qmf_h2g(hcoeff, hrange);
    A0 = A;

    for j=1:size(D,1)
        A0 = wavelet.OrthoWvl.idwt2(A0, D{j,1}, D{j,2}, D{j,3}, hcoeff, hrange, gcoeff);
    end

end