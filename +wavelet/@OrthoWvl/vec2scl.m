function [A, D] = vec2scl(obj, W)
% [A, D] = vec2scl(obj, W)
% Convert a wavelet coefficient vector to [A,D] form, with A, D the approximation and
% detail coefficients.
%
    A.coeff = reshape(W(1:obj.wptr(1)-1), size(obj.wmask.A));
    A.rangex = obj.wrangex.A;
    A.rangey = obj.wrangey.A;

    j = 1; n = 1;
    while n<length(obj.wptr) && obj.wptr(n)-1 < length(W)
        for k=1:3
            if n<length(obj.wptr)
                toto = W(obj.wptr(n):obj.wptr(n+1)-1);
            else
                toto = W(obj.wptr(n):end);
            end
            D{j,k}.coeff = reshape(toto, size(obj.wmask.D{j,k}));
            D{j,k}.rangex = obj.wrangex.D{j,k};
            D{j,k}.rangey = obj.wrangey.D{j,k};
            n = n+1;
        end
        j = j+1;
    end
end
