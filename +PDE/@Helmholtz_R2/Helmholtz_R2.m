classdef Helmholtz_R2 < PDE.Small_Inclusions
% Class for modeling the Helmholtz equation with small inclusion(s).
%     
% Reference [1]: Shape identification and classification in echo-location.
    
    properties(SetAccess = protected)
        pmeb_bg % background magnetic parameter (permeability), or epsilon_0
        pmtt_bg % background electric parameter (permittivity), or mu_0
        
        wavenb_bg % background wave number k0
        wavenb % wave number k
    end
    
    methods
        function obj = Helmholtz_R2(D, cfg, pmeb_bg, pmtt_bg)
            obj = obj@PDE.Small_Inclusions(D, cfg);
            obj.pmeb_bg = pmeb_bg;
            obj.pmtt_bg = pmtt_bg;
        end     
        
        function [F, F_bg, SX, SY, mask] = calculate_field(obj, freq, xs, width, N)
        % Calculate the background potential field and the potential field due to
        % the inclusion.
        % Inputs:
        % s: index of the source
        % xlim: interval in x-axis of rectangular region
        % ylim: interval in y-axis of rectangular region
        % dh: sampling step
        % Outputs:
        % F: potential field u
        % F_bg: potential field U
        % SX, SY: coordinates in x and y-axis of the rectangular region
                        
            if obj.nbIncls > 1
                error('Multiple inclusions are not supported in current version!');
            else
                D = obj.D{1};
                [SX, SY, mask] = D.interior_mesh(max(width, obj.aradius*2*1.5), N);
                
                [F, F_bg] = PDE.Helmholtz_R2.evaluate_field(D, freq, xs, obj.pmeb_bg, obj.pmtt_bg, ...
                                                            SX, SY, mask);
            end
        end
        
        out = data_simulation(obj, freq)

        plot_field(obj, F, F_bg, SX, SY, nbLine, varargin)        
    
        out = reconstruct_SCT(obj, MSR, ord, maxiter, tol)
        out = reconstruct_SCT_analytic(obj, MSR, ord)

        function val = get.wavenb_bg(obj)
            val = obj.freq * sqrt(obj.pmeb_bg*obj.pmtt_bg);
        end
        
        function val = get.wavenb(obj)
            val = obj.D{1}.wavenb(obj.freq);
        end

        function A = MSR2FFP(obj, V, sfreq, mask)            
        % Convert MSR to far field pattern
            [M,N] = size(V{1});
            L = length(V);
            
            A = zeros(M,N,L);

            for n=1:L
                k0 = sfreq(n) * sqrt(obj.pmeb_bg*obj.pmtt_bg);
                cst = sqrt(8*pi*k0*obj.cfg.radius_rcv);
                
                if nargin<4 || isempty(mask)
                    A(:,:,n) = V{n} * cst;
                else
                    A(:,:,n) = V{n} .* mask * cst;
                end
            end
        end
    end
    
    methods(Static)
        function W = post_processing(W0, z0, k0)
        % Post-processing operator for the reconstruction of SCT coefficients.
        % Remark: this operator can induce important error if z0 is not close to the origin due
        % to the truncation (in both W0 and the cylindrical waves).
            
            ord = max(100, size(W0, 1)); % truncation length for the filter
            G=zeros(2*ord+1,1);
            
            for n=-ord:ord
                G(n+ord+1)=tools.Helmholtz.cylind_wave(k0, n, -z0);
            end
            W = conv2(G, 1, W0, 'same');

            % toto = conv2(G, 1, W0); 
            % W = wkeep(toto, size(W0)); % keep the central part
        end
        
        function M = system_matrix_block_A11(freq, D, type, step)
            k = D.kvalH(freq);
            SH = ops.SingleLayer_H(k, D, type, step);
            M = SH.stiffmat;           
        end
        
        function M = system_matrix_block_A12(freq, D, type, step, pmeb_bg, pmtt_bg)
            k0 = freq*sqrt(pmtt_bg * pmeb_bg);
            SH = ops.SingleLayer_H(k0, D, type, step);
            M = -1 * SH.stiffmat; 
        end
        
        function M = system_matrix_block_A21(freq, D, type, step)
            k = D.kvalH(freq);
            Id = ops.Ident(D, type, step);
            KsH = ops.Kstar_H(k, D, type, step);
            M = (-1/2*Id.stiffmat + KsH.stiffmat)/D.pmeb;
        end
        
        function M = system_matrix_block_A22(freq, D, type, step, pmeb_bg, pmtt_bg)
            k0 = freq*sqrt(pmeb_bg*pmtt_bg);
            Id = ops.Ident(D, type, step);
            KsH = ops.Kstar_H(k0, D, type, step);
            M = - (1/2*Id.stiffmat + KsH.stiffmat)/pmeb_bg;
        end
        
        function [b] = source_vector(D, freq, Xi, pmeb_bg, pmtt_bg)
        % Compute the vector of right hand side in eq. 5.11 of [1], namely, U and \partialU\partial\nu,
        % where U is the plane wave: e^{k*i*<\xi,x>}
            
            k0 = freq*sqrt(pmtt_bg*pmeb_bg);
            s1 = tools.Helmholtz.plane_wave(k0, D.points, Xi);
            s2 = tools.Helmholtz.plane_wave_Dn(k0, D, Xi)/pmeb_bg;

            b = [s1; s2]; % Concatenuation gives a 2N-by-M, with N=D.nbPoints, M=number of
                          % source directions
        end
        
        function [b] = source_vector_cw(D, freq, m, pmeb_bg, pmtt_bg)
        % Compute the vector of right hand side in eq. 5.11 of [1], namely, U and \partialU\partial\nu,
        % where U is the cylindrical wave: Jm(k_0|x|)e^{im theta_x}

            k0 = freq*sqrt(pmtt_bg*pmeb_bg);
            s1 = tools.Helmholtz.cylind_wave(k0, m, D.points);
            s2 = tools.Helmholtz.cylind_wave_Dn(k0, m, D)/pmeb_bg;

            b = [s1; s2];
        end
        
        function [sol] = solve_forward(D, freq, Xi, pmeb_bg, pmtt_bg)
        % Solve the forward problem 5.11 of [1] using plane wave as sources U.
            
            if length(freq)>1
                error('Frequency must be a scalar.');
            end

            % Construct the blocks of the system matrix
            matrix_A = PDE.Helmholtz_R2.system_matrix_block_A11(freq, D, 'P0', 1);
            matrix_B = PDE.Helmholtz_R2.system_matrix_block_A12(freq, D, 'P0', 1, pmeb_bg, pmtt_bg);
            matrix_C = PDE.Helmholtz_R2.system_matrix_block_A21(freq, D, 'P0', 1);
            matrix_D = PDE.Helmholtz_R2.system_matrix_block_A22(freq, D, 'P0', 1, pmeb_bg, pmtt_bg);

            b_source = PDE.Helmholtz_R2.source_vector(D, freq, Xi, pmeb_bg, pmtt_bg);

            matrix_BEM = [[matrix_A matrix_B] ; [matrix_C matrix_D]] ;
            % cond(matrix_BEM)
            sol = matrix_BEM\b_source;
        end

        % [sol] = solve_forward(D, freq, Xi, pmeb_bg, pmtt_bg)
        % % function sol = solve_forward_fast(matrix_BEM, b)
        % %     sol = matrix_BEM\b;
        % % end
        
        function [sol] = solve_forward_cw(D, freq, m, pmeb_bg, pmtt_bg)
        % Solve the forward problem 5.11 of [1] using cylindrical wave
        % Um(x)=Jm(k_0|x|)e^{imtheta_x} as source.

            % Construct the blocks of the system matrix
            matrix_A = PDE.Helmholtz_R2.system_matrix_block_A11(freq, D, 'P0', 1);
            matrix_B = PDE.Helmholtz_R2.system_matrix_block_A12(freq, D, 'P0', 1, pmeb_bg, pmtt_bg);
            matrix_C = PDE.Helmholtz_R2.system_matrix_block_A21(freq, D, 'P0', 1);
            matrix_D = PDE.Helmholtz_R2.system_matrix_block_A22(freq, D, 'P0', 1, pmeb_bg, pmtt_bg);

            b = PDE.Helmholtz_R2.source_vector_cw(D, freq, m, pmeb_bg, pmtt_bg);

            matrix_BEM = [[matrix_A matrix_B] ; [matrix_C matrix_D]] ;
            sol = matrix_BEM\b;
        end
                
        function [F, F_bg] = evaluate_field(D, freq, xs, pmeb_bg, pmtt_bg, SX, SY, mask)
        % Evaluate the solution field u and U for a given source xs at the points of coordinates given by
        % SX and SY, which should not contain the interior of D.
            
            SXi = SX(find(mask)); SYi = SY(find(mask));
            Zi = [SXi(:) SYi(:)]';
            SXo = SX(find(1-mask)); SYo = SY(find(1-mask));
            Zo = [SXo(:) SYo(:)]';
            
            % Resolution of the system 5.11
            sol = PDE.Helmholtz_R2.solve_forward(D, freq, xs, pmeb_bg, pmtt_bg);
            vphi = sol(1:size(sol,1)/2, :);
            vpsi = sol(size(sol,1)/2+1:end, :);

            k0 = sqrt(pmtt_bg*pmeb_bg)*freq;
            umU = ops.SingleLayer_H.eval(k0, D, vpsi, Zo);
            U = tools.Helmholtz.plane_wave(k0, Zo, xs); 
            uo = umU(:)+U(:);
            
            ks = sqrt(D.pmtt*D.pmeb)*freq;
            ui =  ops.SingleLayer_H.eval(ks, D, vphi, Zi);
            
            F=zeros(size(mask));
            F(find(mask)) = ui;
            F(find(1-mask)) = uo;
            F_bg = tools.Helmholtz.plane_wave(k0, [SX(:) SY(:)]', xs); 
            F_bg = reshape(F_bg, size(mask));
        end
        
        function out = add_white_noise(data, nlvl)
        % Add white noise to simulated data, see the function data_simulation
        % nlvl: noise level
            
            out = data;
            mode = 1; % treat each source independently
            rowmajor = 1;
            
            % add real white noise
            for n=1:length(data.MSR)
                out.MSR_noisy{n} = tools.add_white_noise(real(out.MSR{n}), nlvl, mode, rowmajor);
            end
            
            % add imaginary white noise
            for n=1:length(data.MSR)
                out.MSR_noisy{n} = out.MSR_noisy{n} + 1i * tools.add_white_noise(imag(out.MSR{n}), ...
                                                                  nlvl, mode, rowmajor);
            end            
        end

        out = make_linop_SCT(cfg, k0, ord)  % construct the linear operator - SCT

        % construct the acquisition matrices of source and receiver
        Ar = make_matrix_Rcv(k0, rcv, z0, ord)
        As = make_matrix_Src(k0, src, z0, ord)
        % [As, Ar] = make_matrix_SR(theta_Xi, k_0, Xr, ord);
    
    end    
end

% function out = data_simulation(obj, freq)
% % Simulation of the MSR matrix.
%     if obj.nbIncls > 1
%         error('NotImplemented: data_simulation for multi-inclusions.')
%     else   
%         if nargin < 2
%             freq = obj.freq;
%         end

%         rcv = obj.cfg.all_rcv;
%         src = obj.cfg.all_src;

%         for n=1:length(obj.freq)
%             out.MSR{n} = PDE.Helmholtz_R2.MSR_simulation(obj.D{1}, freq(n), obj.pmtt_bg, ...
%                                                          obj.pmeb_bg, rcv, src);
%         end
%     end
% end





%         function [sol] = solve_Helm_sys_pwv(D, freq, X, theta_Xi, pmeb_bg, epsilon_0)
%         %Solve the system of [\phi,\psi] in the Helmholtz equation, for each value of frequency, 
%         %we have corresponding pair of solution, we took U is the plane wave function%

%             type1 = 'P0'; step = 1; 

% %             Psi_int = D.sigma(:); 
% %             Phi_int = Psi_int; 

%             % Construct the blocks of the system matrix
%             matrix_A = PDE.Helmholtz_R2.system_matrix_block_A11(freq, D, type1, step);
%             matrix_B = PDE.Helmholtz_R2.system_matrix_block_A12(freq, D, type1, step, pmeb_bg, epsilon_0);
%             matrix_C = PDE.Helmholtz_R2.system_matrix_block_A21(freq, D, type1, step);
%             matrix_D = PDE.Helmholtz_R2.system_matrix_block_A22(freq, D, type1, step, pmeb_bg, epsilon_0);


%             b_source = PDE.Helmholtz_R2.b_vector(D, freq, X, theta_Xi, pmeb_bg, epsilon_0);

% %             matrix_BEM = [[matrix_A matrix_B] ; [matrix_C matrix_D]; 
% %                           [reshape(Psi_int, 1, []), zeros(1,size(matrix_A,2))]; [zeros(1,size(matrix_A,2)), ...
% %                                 reshape(Phi_int,1,[])]] ;
% %             
% %             sol = matrix_BEM\[b_source; zeros(1, size(b_source,2)); zeros(1, size(b_source,2))] ;

%             matrix_BEM = [[matrix_A matrix_B] ; [matrix_C matrix_D]] ;
%             sol = matrix_BEM\b_source ;
%         end        
