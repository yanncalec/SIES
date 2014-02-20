function mask = mask_onescl(obj, j, k)
% Make a 0-1 mask for extraction of wavelet coefficients of a given scale and
% direction. This mask is different to the mask obj.wmask which selects the active
% coefficients and applied only on one scale and direction, and it applies on the
% whole coefficient vector.
    
    mask = zeros(size(obj.wmask.all));
    if j==0
        mask(1:obj.wptr(1)-1) = 1;
    else
        nn = 3*(j-1)+k;
        if nn==length(obj.wptr)
            mask(obj.wptr(nn):end) = 1;
        else
            mask(obj.wptr(nn):obj.wptr(nn+1)-1) = 1;
        end
    end            
end
