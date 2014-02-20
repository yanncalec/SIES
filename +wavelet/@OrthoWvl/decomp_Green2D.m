function [Wa, Xim] = decomp_Green2D(obj, Z)
% Calculate wavelet coefficients of the Green function of Laplacian.  The function
% f(X)=1/2pi*log|X-Z| is evaluated on the AROI grid and the wavelet coefficients are
% computed using fast transform algorithms. Only those active, (i.e. which contribute
% to the ROI) are kept.
%
% INPUT:
% Z: singularity of the Green function, a 2d vector
% OUTPUTS:
% Wa: coefficients of active wavelets
% Xim: value of the Green function (times the constant 2^Jnum) evaluated on AROI grid, an image

    
    % % Verify the singularity is outside the AROI
    % if (Z(1)>=obj.AROI.suppx(1) && Z(1)<=obj.AROI.suppx(2)) && (Z(2)>=obj.AROI.suppy(1) ...
    %                                                       && Z(2)<=obj.AROI.suppy(2))
    %     warning('The singularity is included inside the augmented ROI!');
    %     fprintf('Z(1)=%f, AROIx=[%f, %f]\n', Z(1),obj.AROI.suppx(1), ...
    %             obj.AROI.suppx(2));
    %     fprintf('Z(2)=%f, AROIy=[%f, %f]\n', Z(2),obj.AROI.suppy(1), ...
    %             obj.AROI.suppy(2));
    % end
    
    % The point value of the function approximates the scaling coefficients
    % (up to a factor) in very fine scales:
    
    % the singularity of Green function is outside the ROI hence removed:
    X = tools.Laplacian.Green2D_trunc(obj.AROI.meshpoints, Z(:)) * 2^obj.Jnum;     
    % X = tools.Laplacian.Green2D(obj.AROI.meshpoints, Z(:)) * 2^obj.Jnum;

    Xim = reshape(X, obj.AROI.dim);
    [~, ~, W] = obj.analysis(Xim); % decomposition from Jnum to Jmax
    Wa = W .* obj.wmask.all; % set inactive wavelets to zero
end
