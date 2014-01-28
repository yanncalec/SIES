function out = reconstruct_PT(obj, SFR, Current_bg, maxiter, tol, symmode)
% Reconstruct PT from the post-processed dipolar expansion SFR matrix.
% The output is a structure including the following fields:
% As, Ar: the left and right side matrix of the forward operator
% PT: reconstructed PT matrix
% res: residual error of |L(PT) - V)|
% Lsymm, Lnsym: linear operators with and without the constraint of
% symmetry

    if nargin < 6
        symmode = 0;
    end
    if nargin < 5
        tol = 1e-3;
    end
    if nargin < 4
        maxiter=100;
    end

    % Construct the linear operator L
    op = PDE.ElectricFish.make_linop_PT(obj.cfg, Current_bg, obj.impd, symmode);

    [X, ~] = lsqr(op.L, SFR(:), tol, maxiter);
    
    PT = reshape(X, 2, 2);
    out.res = norm(SFR(:) - op.L(PT, 'notransp'));
    out.op = op; 

    if symmode
        PT = PT + PT.';
    end
    out.PT = PT;
end    

