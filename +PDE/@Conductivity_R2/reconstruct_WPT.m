function out = reconstruct_WPT(MSR, Lop, mask, maxiter, tol)
% Simple least square reconstruction of the wavelet coefficients
    import tools.linsys.* PDE.Conductivity_R2

    if nargin < 5
        tol = 1e-5;
    end
    if nargin < 4
        maxiter = 10^4;
    end

    if nargin<3 || isempty(mask)
        Lopm = Lop;
    else
        Lopm = @(X, tflag)tools.linsys.linop_mask(X, Lop, mask(:), tflag);
    end

    out.X = lsqr(Lopm, MSR(:), tol, maxiter);        
end
