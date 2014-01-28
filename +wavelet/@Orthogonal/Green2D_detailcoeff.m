function D = Green2D_detailcoeff(obj, Z, nbScl)
    
    D = cell(nbScl, 3);
    
    for j=1:nbScl
        J = obj.Jmax-j+1;

        D{j,1}.coeff = wavelet.OrthoWvl.Green2D_detailcoeff_mex(obj.Tphi, obj.phisupp, obj.Tpsi, ...
                                                          obj.psisupp, J, obj.DetailSpace{j,1}.rangex, ...
                                                          obj.DetailSpace{j,1}.rangey, Z);
        D{j,2}.coeff = wavelet.OrthoWvl.Green2D_detailcoeff_mex(obj.Tpsi, obj.psisupp, obj.Tphi, ...
                                                          obj.phisupp, J, obj.DetailSpace{j,2}.rangex, ...
                                                          obj.DetailSpace{j,2}.rangey, Z);
        D{j,3}.coeff = wavelet.OrthoWvl.Green2D_detailcoeff_mex(obj.Tpsi, obj.psisupp, obj.Tpsi, ...
                                                          obj.psisupp, J, obj.DetailSpace{j,3}.rangex, ...
                                                          obj.DetailSpace{j,3}.rangey, Z);
    end        
end
