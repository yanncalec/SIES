function out = reconstruct_CGPT(obj, MSR, Current, ord, maxiter, tol, symmode)
% Reconstruct the contracted GPT from data using iterative method.
%
% The output is a structure including the following fields:
% As, Ar: the left and right side matrix such that
%         MSR = L(CGPT) = As * CGPT * transpose(Ar).
%         In case of equispaced circular full view setting, As and Ar are
%         written as Cs * Ds and Cr * Dr (see the function make_matrix_CD). In
%         general case these are linear operators
% CGPT: reconstructed CGPT matrix, an object of GPT.CGPT class
% res: residual error of |L(CGPT) - MSR)|

if nargin < 7
	symmode = 0;
end
if nargin < 6
	tol = 1e-5;
end
if nargin < 5
	maxiter = 10^5;
end

% Construct the linear operator L
op = PDE.ElectricFish.make_linop_CGPT(obj.cfg, Current, obj.impd, ord, symmode);

[X, ~] = lsqr(op.L, MSR(:), tol, maxiter); % LSQR method

CGPT = reshape(X, 2*ord, 2*ord);
out.res = norm(MSR(:) - op.L(CGPT, 'notransp'));
out.op = op;

if symmode
	CGPT = CGPT+CGPT.';
end
out.CGPT = CGPT;
end
