function [Xr, err] = apprx_ROI(obj, W, nbScl, Xim)
% Resynthesis of a function on the ROI by using coefficients of coarse scales.
% INPUTS:
% W: wavelet coefficients
% nbScl: number of detail spaces used for synthesis, 0 for approximation space only. By default
% all scales are used.
% Xim: original function on the AROI, optional
% OUTPUTS:
% Xr: resynthesized function on the ROI
% err: approximation error

    if nargin >= 3
        W = obj.extract_coeffs_lowscls(W, nbScl); % approximation coefficients
    end

    [A, D] = obj.vec2scl(W);
    Xr = obj.synthesis_ROI(A,D);

    if nargin == 4
        X0 = obj.ROI.resize(Xim, obj.AROI, 'crop'); % Restrict on the ROI
        err = norm(Xr - X0, 'fro');
    else 
        err = [];
    end
end
