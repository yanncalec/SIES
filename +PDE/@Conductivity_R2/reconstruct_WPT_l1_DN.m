function out = reconstruct_WPT_l1_DN(MSR, Lop, W, M, rho, maxiter, tol)
% Reconstruct the WPT coefficients by solving a l1 minimization (denoising) problem:
% min_X 1/2*||L(M.*X)-V||^2 + rho * ||W.*X||_1
% INPUTS:
% MSR: the data matrix V
% Lop: linear operator
% W: weight of the L1 norm, eg, the norm of column vectors of the forward operator
% M: 0-1 mask of the forward operator
% rho: threshold parameter, e.g., rho <= sigma_noise * sqrt(2*log(#M))
% maxiter: maximal number of iterations for l1 algorithm
% tol: tolerance of the l1 algorithm

if nargin < 7
	tol = 1e-5;
end
if nargin < 6
	maxiter = 10^4;
end

if ~isempty(M)
	Lopm = @(X, tflag)tools.linsys.linop_mask(X, Lop, M(:), tflag);
else
	Lopm = Lop;
end

% L1-Fista
out = tools.solver.l1fista(Lopm, MSR(:), W, rho, tol, maxiter, 1000); % print message every
% 1000 iterations


% % YALL1
% opts.tol = tol;
% opts.maxiter = maxiter;
% opts.rho = rho;
% opts.nonortho = 1;
% opts.weights = W(:);
% opts.print = 2;
% [out.X, out.msg] = tools.solver.yall1(Lopm, MSR(:), opts);

end