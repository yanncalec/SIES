function out = reconstruct_SCT(obj, MSR, freq, ord, maxiter, tol)

if nargin < 5
	tol = 1e-5;
end
if nargin < 4
	maxiter = 10^5;
end

% ord = min(ord, min(floor((obj.cfg.Ns-1)/2), floor((obj.cfg.Nr-1)/2))) ; % maximum order

op = PDE.Helmholtz_R2.make_linop_SCT(obj.cfg, obj.wavenb_bg(freq), ord); % Construct the linear operator L
L = op.L;

X = lsqr(L, MSR(:), tol, maxiter); % LSQR method

out.SCT = reshape(X, 2*ord+1, 2*ord+1);
out.res = norm(MSR(:) - L(out.SCT, 'notransp'));
out.op = op;

% out.SCT0 = SCT0;
% out.SCT = PDE.Helmholtz_R2.post_processing(SCT0, obj.cfg.center, obj.wavenb_bg);

end

