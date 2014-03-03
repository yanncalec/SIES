function WPT = theoretical_WPT_std(D, tvec, normal, avec, sigma, kappa, Z0, Jmax, ...
                                   ApprxSpace, DetailSpace, Tpsi, Tphi, Gpsi, ...
                                   Gphi)
    % Compute the WPT matrix of standard representation.
    %
    % INPUTS:
    % OUTPUTS:
    % WPT: an object of WPT_std class
    %
    % Convention:
    % The row index m in the WPT coefficient matrix W(m,n) corresponds to the
    % "source", ie, the second function psi in 
    % W(m,n) = \int_{\p D} psi_n (\lambda I-K_D^*)^-1 [d psi_m dn] ds(x)
    %

    nbPoints = size(D,2);
        
    % Threshold of truncation. Coefficients smaller than this value are ignored.
    trunc_epsilon = 1e-10;
    
    nbDir = 3;

    % Number of detail spaces
    nbScl = size(DetailSpace, 1); 
    Jmin = Jmax - nbScl + 1;

    % The inverse operator of (lambda I - K_D*)
    lambda = (kappa+1)/2/(kappa-1);
    Ks = ops.Kstar.make_kernel_matrix(D, tvec, normal, avec, sigma);
    A = lambda*eye(nbPoints, nbPoints) - Ks;
    %[L,U] = lu(A); % LU decomposition of A
    Ainv = inv(A);

    clear A Ks;

    %% Evaluation of wavelets on the boundary: Approximation space
    N =  ApprxSpace.pmeshpoints;
    [Vr, dX, dY] = wavelet.OrthoWvl.evaluation(D, Jmax, N, Tphi, Tphi, Gphi, Gphi);
    dVs = diag(normal(1,:)) * dX + diag(normal(2,:)) * dY;

    clear dX dY;

    %% Computation of wavelet coefficients: Apprx-Apprx
    AA = tools.full2sparse(transpose(Ainv*dVs) * (diag(sigma) * Vr), trunc_epsilon);

    %% Detail coefficients exist if nbScl>0
    if nbScl>0
        %% Evaluation of wavelets on the boundary: Detail spaces
        Wr = cell(nbDir, nbScl);
        dWs = cell(nbDir, nbScl);
        
        for l=1:nbDir
            for J=Jmin:Jmax 
                N = DetailSpace{Jmax-J+1, l}.pmeshpoints;
                % Wavelets for receptor of this scale
                if l==1
                    [val, dX, dY] = wavelet.OrthoWvl.evaluation(D, J, N, Tphi, Tpsi, Gphi, Gpsi);
                elseif l==2
                    [val, dX, dY] = wavelet.OrthoWvl.evaluation(D, J, N, Tpsi, Tphi, Gpsi, Gphi);
                else                
                    [val, dX, dY] = wavelet.OrthoWvl.evaluation(D, J, N, Tpsi, Tpsi, Gpsi, Gpsi);
                end                
                j = J-Jmin+1;
                Wr{l,j} = sparse(val);
                dWs{l,j} = sparse(diag(normal(1,:)) * dX + diag(normal(2,:)) * dY);
            end
        end

        clear dX dY;

        %% Computation of wavelet coefficients
        DD = cell(nbDir);
        DA = cell(nbDir,1);
        AD = cell(1,nbDir);
        
        % Detail-Detail (or D)
        %
        % Loops on the direction
        for lr=1:nbDir
            for lc=1:nbDir
                DD{lr,lc} = cell(nbScl);
                % Loops on the scale
                for jr=1:nbScl
                    for jc=1:nbScl                                            
                        DD{lr,lc}{jr,jc} = tools.full2sparse(transpose(Ainv*dWs{lr,jr})* ...
                                                             diag(sigma)*Wr{lc,jc}, trunc_epsilon);                        
                        
                    end
                end
            end
        end

        % Detail-Apprx (or B)
        %
        % Loops on the direction
        for lr=1:nbDir
            DA{lr} = cell(nbScl,1);
            % Loops on the scale    
            for jr=1:nbScl            
                DA{lr}{jr} = tools.full2sparse(transpose(Ainv*dWs{lr,jr}) * (diag(sigma) * Vr), ...
                                               trunc_epsilon);
                
            end
        end

        % Apprx-Detail (or C)
        %
        % Loops on the direction
        for lc=1:nbDir
            AD{lc} = cell(1,nbScl);
            % Loops on the scale
            for jc=1:nbScl
                AD{lc}{jc} = tools.full2sparse(transpose(Ainv*dVs) * (diag(sigma) * Wr{lc,jc}), trunc_epsilon);
            end
        end
        
        WPT = wavelet.WPT_std(DD, DA, AD, AA);    
    else
        WPT = wavelet.WPT_std({}, {}, {}, AA);    
    end
end