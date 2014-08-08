function out = reconstruct_CGPT(obj, MSR, ord, maxiter, tol, symmode, method)
% Reconstruct the contracted GPT from data using iterative method.
%
% Inputs:
% MSR: the MSR data matrix, if it is time-dependent or multifrequency, MSR
% need to be a cell.
% ord: maximum order of reconstruction
% maxiter, tol: see lsqr method, not needed for pinv method
% symmmode: if true a symmetric constraint on the solutionwill be incoorperated
% method: 'pinv' for formal inversion by pseu-inverse or'lsqr' by iteration. 
%
% Output:
% A structure including the following fields:
% As, Ar: the left and right side matrix such that
%         MSR = L(CGPT) = As * CGPT * transpose(Ar).
%         In case of equispaced circular full view setting, As and Ar are
%         written as Cs * Ds and Cr * Dr (see the function make_matrix_CD). In
%         general case these are linear operators
% CGPT: reconstructed CGPT matrix, an object of GPT.CGPT class
% res: residual error of |L(CGPT) - MSR)|

if ~iscell(MSR) % Convert to a cell, for compability with PulseImaging_R2 class
    MSR = {MSR};
end

if nargin < 7
    method = 'pinv';
end

if nargin < 6
    symmode = 0;
end
if nargin < 5
    tol = 1e-5;
end
if nargin < 4
    maxiter = 10^5;
end

op = PDE.Conductivity_R2.make_linop_CGPT(obj.cfg, ord, symmode); % Construct the linear operator L
out.op = op;

if strcmp(method, 'lsqr')
    for t=1:length(MSR)
        toto = reshape(MSR{t}, [], 1);
        X = lsqr(op.L, toto, tol, maxiter); % LSQR method
        
        CGPT = reshape(X, 2*ord, 2*ord);
        out.res{t} = norm(toto - op.L(CGPT, 'notransp'));
        
        if symmode
            CGPT = CGPT+CGPT.';
        end
        out.CGPT{t} = CGPT;
    end
elseif strcmp(method, 'pinv')
    for t=1:length(MSR)
        CGPT = pinv(op.As)*MSR{t}*pinv(op.Ar');
        out.res{t} = norm(MSR{t} - op.As*CGPT*op.Ar', 'fro');        
        out.CGPT{t} = CGPT;
    end
    
end

