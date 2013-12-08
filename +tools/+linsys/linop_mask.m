function Y=linop_mask(X, op, mask, tflag)
% Modify a linear operator by applying a 0-1 mask to it
% INPUTS:
% X: input vector
% op: matlab linear operator (see cgs, lsqr for example)
% mask: 0-1 mask

    if strcmp(tflag,'notransp')
        Y = op(X.*mask, 'notransp');
    elseif strcmp(tflag,'transp')
        Y = mask.*op(X, 'transp');
    end
end