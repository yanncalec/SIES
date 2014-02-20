function mask = mask_active_lowscls(obj, nbScl)            
% Make a 0-1 mask for extraction of wavelet coefficients up to a given scale

    mask = obj.mask_active_onescl(0,0);            
    for j=1:nbScl
        mask = obj.mask_active_onescl(j,1) + obj.mask_active_onescl(j,2) + ...
               obj.mask_active_onescl(j,3) + mask;
    end
end

