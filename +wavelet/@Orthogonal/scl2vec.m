function W = scl2vec(A, D, nbScl)
% W = scl2vec(A, D, nbScl)
% Keep the first nbScl scales in the wavelet coeffcients A, D returned by wavedec2 and
% concatenate coefficients into a vector form
%
    if nargin < 3
        nbScl = size(D,1);
    end

    W = A.coeff(:); 

    for j=1:nbScl
        W = [W; D{j,1}.coeff(:); D{j,2}.coeff(:); D{j,3}.coeff(:)];
    end            
end
