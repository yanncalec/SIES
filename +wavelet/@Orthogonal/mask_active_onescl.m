function mask = mask_active_onescl(obj, j, k)
    mask = zeros(size(obj.wmask.all));
    if j==0
        mask(1:obj.wptr(1)-1) = obj.wmask.A(:);
    else
        nn = 3*(j-1)+k;
        if nn==length(obj.wptr)
            mask(obj.wptr(nn):end) = obj.wmask.D{j,k}(:);
        else
            mask(obj.wptr(nn):obj.wptr(nn+1)-1) = obj.wmask.D{j,k}(:);
        end
    end            
end
