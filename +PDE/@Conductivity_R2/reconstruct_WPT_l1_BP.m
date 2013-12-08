function out = reconstruct_WPT_l1_BP(MSR, Lop, W, M, delta, maxiter, tol)
% Reconstruct the WPT coefficients by solving a l1 minimization (basis-pursuit) problem:
% min_X ||W.*X||_1 st ||L(M.*X)-V|| <= delta
% INPUTS:
% MSR: the data matrix V
% Lop: linear operator
% W: weight of the L1 norm, eg, the norm of column vectors of the forward operator
% M: 0-1 mask of the forward operator
% delta: noise level
% maxiter: maximal number of iterations for l1 algorithm
% tol: tolerance of the l1 algorithm

    if nargin < 7
        tol = 1e-5;
    end
    if nargin < 6
        maxiter = 10^4;
    end
    if nargin < 5
        delta = 0;
    end

    if delta > 0
        error('Basis pursuit denoising not implemented!');
    end
    
    if ~isempty(M)
        Lopm.times = @(X)tools.linsys.linop_mask(X, Lop, M(:), 'notransp');
        Lopm.trans = @(X)tools.linsys.linop_mask(X, Lop, M(:), 'transp');
    else
        Lopm.times = @(X)Lop(X,'notransp');
        Lopm.trans = @(X)Lop(X,'transp');
    end
    
    % YAll1
    opts.tol = tol;
    opts.maxiter = maxiter;
    opts.delta = delta;
    opts.nonortho = 1;
    opts.weights = W(:);
    opts.print = 2;

    [out.X, out.msg] = tools.solver.yall1(Lopm, MSR(:), opts);
    
end
