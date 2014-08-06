classdef Conductivity_R2 < PDE.Small_Inclusions
% Class for the conductivity problem in the free space
% u solves:
% div(rho_D*grad(u)) = 0, in R^2 \ x_s
% u-G = O(1/x) as |x|->oo
% where rho_D = 1 + \sum_n (k_n - 1)indic(D_n), with k_n = cnd_n + 1i * freq * pmtt_n,
% and x_s is the position of the source, and G(x) = 1/2/pi*ln|x|.
%
% This class can only simulate data with the mono-Dirac source. For the multiple-Dirac (dipole source) case, one should modify
% the function compute_dGdn, see the class |PulseImaging_R2| for more details.

    properties(Access=protected)
        % Relative to the linear system A[phi] = b.
        % A is a block matrix: (lambda I - K_{D_n}^*) on the diagonal, dSdn on the off diagonal
        % b can be the normal derivative of either the polynomials (for computation of GPTs), or
        % the Green's function (for simulation of data)
        % phi is the solution and it is used in either the computation of GPTs or the simulation
        % of MSR: u-G = \sum_l dSdn[phi_l]
        %
        % In the matrix A only the diagonal terms depend on the frequency (via lambda I), the
        % other frequency-independent coefficients are stored in KsdS.
        % We store also dGdn which is frequency-independent.
        KsdS=[] % system block matrix
        dGdn=[] % the normal derivative of the Green's function, also the right hand vector in A[phi]=b
    end
    
    properties(SetAccess=protected)
        cnd % conductivity constant corresponding to each inclusion, an array
        pmtt % permittivity constant corresponding to each inclusion, an array
    end
    
    methods(Access=protected) % Auxiliary functions
        function val = compute_dGdn(obj, s)
        % Construct the right hand vector of the given source.
        % If the source is not given, compute for all sources
            
            if nargin<2
                src = obj.cfg.all_src();
            else
                src = obj.cfg.src(s);
            end
            
            nbPoints = obj.D{1}.nbPoints;
            val = zeros(nbPoints*obj.nbIncls, size(src,2)); 
            idx=0;
            for i=1:obj.nbIncls
                toto = tools.Laplacian.Green2D_Dn(src, obj.D{i}.points, obj.D{i}.normal);
                val(idx+1:idx+nbPoints, :) = toto';
                idx = idx+nbPoints; 
            end    
        end
        
        function Phi = compute_phi(obj, freq, s)
        % Compute the boundary function phi given the frequency and the source index
        % The function phi_s (s is the index of source) is the solution to 
        % A[phi_s] = b_s
        % A is independent of s, we solve for multiple s by
        % [phi_1, phi_2, ...] = A\[b_1, b_2, ...]
            
            nbPoints = obj.D{1}.nbPoints;
            lambda = asymp.CGPT.lambda(obj.cnd, obj.pmtt, freq);

            % Generate the system matrix
            Amat = asymp.CGPT.make_system_matrix_fast(obj.KsdS, lambda);

            if nargin<3
                dGdn = obj.dGdn;                
            else
                dGdn = obj.compute_dGdn(s);
            end

            phi = Amat\dGdn; % function phi of all sources and inclusions
            Phi = cell(obj.nbIncls, 1);
            
            idx=0;                        
            for i=1:obj.nbIncls
                Phi{i} = phi(idx+1:idx+nbPoints,:); % function phi of the i-th inclusion
                idx = idx+nbPoints;
            end
        end
    end

    methods
        function obj = Conductivity_R2(D, cnd, pmtt, cfg)
            obj = obj@PDE.Small_Inclusions(D, cfg);

            if length(cnd)<obj.nbIncls || length(pmtt)<obj.nbIncls
                error('The value of conductivity and permittivity must be specified for each inclusion!');
            end
            
            for n=1:obj.nbIncls
                if cnd(n)==1 || cnd(n)<0
                    error('The conductivity constant must be positive and different from 1!');
                end
                
                if pmtt(n)<0
                    error('The permittivity constant must be positive!');
                end                
            end
            
            obj.cnd = cnd; obj.pmtt = pmtt;
            
            obj.KsdS = asymp.CGPT.make_block_matrix(obj.D); 
            obj.dGdn = obj.compute_dGdn();
        end
        
        out = data_simulation(obj, freq)

        % Calculate and plot the potential fields 
        [F, F_bg, Sx, Sy, mask] = calculate_field(obj, freq, s, z0, width, N)        
        plot_field(obj, s, F, F_bg, Sx, Sy, nbLine, varargin)
        
        % Reconstruction of contracted GPT matrix
        out = reconstruct_CGPT(obj, MSR, ord, maxiter, tol, symmode, method) 
        out = reconstruct_CGPT_analytic(obj, MSR, ord)
        
        % Reconstruction of wavelet polarization tensors
        out = reconstruct_WPT_l1_DN(MSR, Lop, W, M, rho, maxiter, tol) % L1 denoising
        out = reconstruct_WPT_l1_BP(MSR, Lop, W, M, delta, maxiter, tol) % L1 basis-pursuit        
        out = reconstruct_WPT(MSR, Lop, mask, maxiter, tol)                
    end
        
    methods(Static) % Utility functions                            
        A = make_matrix_A(Xs, Z, order)
        out = make_linop_CGPT(cfg, ord, symmode)  % construct the linear operator - CGPT
        
        A = make_matrix_WPT(Xs, Z, WF, nbScl, stdmode)
        out = make_linop_WPT(cfg, WF, nbScl, stdmode) % construct the linear operator - WPT
                
        function out = add_white_noise(data, nlvl)
        % Add white noise to simulated data.
        % nlvl: noise level
            
            out = data;
            mode = 1; % treat each source independently
            rowmajor = 1; % procede row-by-row
            
            out.sigma = zeros(1, length(data.MSR));
            
            for f=1:length(data.MSR) % f can be the index of frequency or time                
                % add real white noise
                [totor, sigmar] = tools.add_white_noise(real(out.MSR{f}), nlvl, mode, rowmajor); 
                
                % add imaginary white noise
                [totoi, sigmai] = tools.add_white_noise(imag(out.MSR{f}), nlvl, mode, rowmajor);
                toto = totor + 1i * totoi; 

                out.sigma(f) = abs(sigmar + 1i * sigmai);
                out.MSR_noisy{f} = toto;
            end
        end
    end
end
