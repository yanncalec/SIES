function val = extract_coeffs(obj, W, j, k)
% Extract from the vector of wavelet coefficients those of a given scale and
% direction, and reshape to a matrix.
% INPUTS:
% W: input vector        
% j: scale index, j=0 extracts the approximation coeffs
% k: k=1,2,3 for j>0, no effect for j=0

    if nargin < 4
        k = 0;
    end

    mask = obj.mask_onescl(j,k);
    val = W(find(mask));

    if j==0
        val = reshape(val, size(obj.wmask.A));
    else
        val = reshape(val, size(obj.wmask.D{j,k}));
    end
end
