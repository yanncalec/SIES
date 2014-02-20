function W = invert_scl(obj, A, D)
% Invert the scale order in wavelet coefficients.
% INPUTS:
% A,D: approximation and detail coefficients. D{end,1}..D{end,3} are the finest
% scale (Jmin), while D{1,1}..D{1,3} are the coarsest scale.
% OUTPUT
% W: the vector of concatenation:
% W = [D{end,1}, D{end-1,1}...D{1,1}, D{end,2}...D{1,2}, D{end,3}...D{1,3}, A]

    W=[];
    for j=size(D,1):-1:1
        for k=1:length(D(j,:))
            W = [W; D{j,k}.coeff(:)];
        end
    end
    W = [W; A.coeff(:)];
end
