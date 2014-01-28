function [Wa, A, D, Xim] = decomp_Green2D_num(obj, Z, nbScl)
% Calculate wavelet coefficients of the Green function of Laplacian.  The wavelet coefficients are
% computed by the inner product. Only those active, (i.e. which contribute to the ROI) are kept.
%
% INPUT:
% Z: singularity of the Green function, a 2d vector
% nbScl: number of decomposition scales
% OUTPUTS:
% Wa: coefficients of active wavelets
% Xim: value of the Green function (times the constant 2^Jnum) evaluated on AROI grid, an image
    
    % For the approximation coefficients, evaluate by quadrature    
    A = obj.Green2D_apprxcoeff(Z);    
    
    D = obj.Green2D_detailcoeff(Z, nbScl);

    Wa = obj.scl2vec(A,D);

    X = tools.Laplacian.Green2D(obj.AROI.meshpoints, Z(:)) * 2^obj.Jnum;
    Xim = reshape(X, obj.AROI.dim);
end
