function mask = mask_lowscls(obj, nbScl)            
% Make a 0-1 mask for extraction of wavelet coefficients up to a given scale

    mask = obj.mask_onescl(0,0);            
    for j=1:nbScl
        mask = obj.mask_onescl(j,1) + obj.mask_onescl(j,2) + obj.mask_onescl(j,3) + ...
               mask;
    end
end
