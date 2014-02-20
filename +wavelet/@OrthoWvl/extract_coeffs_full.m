function val = extract_coeffs_full(obj, W, j, k)
% Extract from the vector of wavelet coefficients those of a given scale and direction.
% This function is similar to extract_coeffs, but the extracted vector has the same
% size as the vector W.
% 
% INPUTS:
% W: input vector        
% j: scale index, j=0 extracts the approximation coeffs
% k: k=1,2,3 for j>0, no effect for j=0

    if nargin < 4
        k = 0;
    end

    mask = obj.mask_onescl(j,k);
    val = W .* mask;
end
