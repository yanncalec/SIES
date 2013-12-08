function out = l1fista(A, Y, W, mu, tol, maxIter, verbose)
% L1 minimization by Fast Iterative Soft-Thresholding Algorithm (FISTA).
% Solve the l1 minimization problem:
% min_X 1/2*|A(X)-Y|^2 + mu * |W.*X|_1
% 
% A : system operator
% Y : data in vector form
% W : weight for L1 norm
% mu : L1 fidelity penalty, the most important argument, should be set according to noise level
% tol : iteration tolerance
% maxIter : maximum iteration steps
% verbose: if >0, print message every verbose iterations.
    
    A0 = A;
    if isa(A0, 'numeric')
        A = @(X, tflag) tools.linsys.matrix_as_linop(A0, X, tflag);
    end
    
    X0 = A(Y, 'transp'); 
    X = zeros(size(X0));
    Z = X; 
    Z0 = X;

    AX = A(X, 'notransp');
    AX0 = AX;
    AZ = AX;

    gradZ = zeros(size(X));
    gradZ0 = zeros(size(X));
    dgradZ = zeros(size(X));

    gradstep = 0;		% BB gradient step
    RdX = 1.;
    RdX0= 1.; % Relative change in X, iteration stopping criteria
    niter = 0;	 % Counter for iteration

    nY = norm(Y, 'fro'); % norm of Y
    tau = 1.;
    tau0 =1.;
    
    % Support of X
    XSupp = zeros(size(X));
    XSupp0 = XSupp; %(X.abs()>0).select(1, ArrayXi::Zero(X.size()));
    Xdiff = numel(X);
    rXdiff=1.;		% Relative size of support change

    % % Main support (nnz-largest coeffs)
    % MSupp = zeros(size(X));	
    % MSupp0 = XSupp0;
    % Mdiff = numel(X);
    % rMdiff=1.;		% Relative size of main support change

    converged = false;

    while niter < maxIter && ~converged
        gradZ = A(AZ-Y, 'transp');

        if (niter == 0) 
            % Do steepest descent at 1st iteration
            % Gradient of AL function      
            AgradZ = A(gradZ, 'notransp');
            gradstep = sum(gradZ.*gradZ) / sum(AgradZ.*AgradZ);
        else
            % Set gradstep through BB formula
            dgradZ = gradZ - gradZ0;
            gradstep = sum(dZ .* dZ) / sum(dZ .* dgradZ);
        end
        
        X = tools.l1_shrink(Z - gradstep*gradZ, gradstep*mu, W);

        AX = A(X, 'notransp');
        dX = X - X0;
        tau = (1+sqrt(1+4*tau*tau))/2;
        Z = X + (tau0 - 1)/tau * dX;
        dZ = Z - Z0;
        
        AZ = (1+(tau0 - 1)/tau)*AX - (tau0-1)/tau * AX0;

        X0 = X;
        Z0 = Z;
        AX0 = AX;
        
        gradZ0 = gradZ;
        tau0 = tau;

        if niter==0
            RdX = 1;
        else
            RdX = norm(dX, 'fro') / norm(X0, 'fro');
        end

        dX0 = dX;

        % if mod(niter,sFreq) == 0
        %     % Support changement of X
        %     XSupp = tools.support_mask(X); 
        %     Xdiff = sum(abs(XSupp0 - XSupp));
        %     XSupp0 = XSupp;
        %     rXdiff = Xdiff*1./numel(X);
        %     converged = converged || (rXdiff<tol);
        % end

        converged = converged || (RdX < tol && RdX0 < tol) || tools.l0norm(X)==0; % Must have two consequtive small values
        RdX0 = RdX;

        % Print convergence information
        if verbose>0 && mod(niter, verbose) == 0	
            res = norm(AX-Y, 'fro');
            nL1 = sum(abs(X));
            nL0 = tools.l0norm(X);

            fprintf(['Iter.: %d RdX=%1.2e  |AX-Y|=%1.5e  L1norm=%1.2e  nonzero per.=' ...
            '%1.2e\n'], niter, RdX, res, nL1, nL0/numel(X));
            
            % fprintf(['Iteration : %d\tRdX = %1.2e\t|AX-Y| = %1.5e\tL1 norm = %1.2e\tnon zero per. ' ...
            %          '= %1.2e\tXdiff=%d\trXdiff=%1.2e\n'], niter, RdX, res, nL1, nL0/numel(X), Xdiff, rXdiff);
            
        end
        niter = niter + 1;
    end

    if verbose>0 && ~converged
        fprintf('L1FISTA terminated without convergence.\n');
    end

    out.X = X;
    out.niter= niter;
    out.res = norm(AX-Y, 'fro');
    out.vobj = 0.5 * res * res + mu * nL1;
end
