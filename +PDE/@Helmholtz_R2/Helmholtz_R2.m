classdef Helmholtz_R2 < PDE.Small_Inclusions
	% Class for modeling the Helmholtz equation with small inclusion(s).
	%
	% Reference: Shape identification and classification in echo-location.
	
	properties(SetAccess = protected)
		pmeb_bg % background magnetic parameter (permeability), or epsilon_0
		pmtt_bg % background electric parameter (permittivity), or mu_0
		
		pmtt % permittivity constant of inclusion
		pmeb % permeability constant of inclusion
	end
	
	methods
		function val = wavenb(obj, freq)
			% the frequency dependant wave number k for the inclusion
			val = tools.Helmholtz.wavenb(freq, obj.pmeb, obj.pmtt);
		end
		
		function val = wavenb_bg(obj, freq)
			% the frequency dependant wave number k0 for the background
			val = tools.Helmholtz.wavenb(freq, obj.pmeb_bg, obj.pmtt_bg);
		end
		
		function obj = Helmholtz_R2(D, pmeb, pmtt, pmeb_bg, pmtt_bg, cfg)
			% Inputs:
			% -D: an (or a list of) instance(s) of C2boundary class
			% -pmeb: value(s) of permeability
			% -pmtt: value(s) of permittivity
			% -pmeb_bg: value of background permeability			
			% -pmtt_bg: value of background permittivity
			% cfg: configuration, eg, an instance of acq.Planewave
			
			if iscell(D) && length(D)>1
				error('Multi-inclusions for Helmholtz equation are not supported in the current version!');
			end
			
			obj = obj@PDE.Small_Inclusions(D, cfg);
			
			if (pmtt==pmtt_bg) || (pmtt<0)
				error('The permittivity constant must be positive and different from the background');
			end
			
			if pmeb==pmeb_bg || pmeb<0
				error('The permeability constant must be positive and different from the background');
			end
			
			%             if length(pmtt)<obj.nbIncls || length(pmeb)<obj.nbIncls
			%                 error('The value of permeability and permittivity must be specified for each inclusion');
			%             end
			
			%             for n=1:obj.nbIncls
			%                 if pmtt(n)==pmtt_bg || pmtt(n)<0
			%                     error('The permittivity constant must be positive and different from the background');
			%                 end
			%
			%                 if pmeb(n)==pmeb_bg || pmeb(n)<0
			%                     error('The permeability constant must be positive and different from the background');
			%                 end
			%             end
			
			obj.pmtt = pmtt;
			obj.pmeb = pmeb;
			
			obj.pmtt_bg = pmtt_bg;
			obj.pmeb_bg = pmeb_bg;
		end
		
		function [F, F_bg, SX, SY, mask] = calculate_field(obj, freq, s, z0, width, N)
			% Calculate the background potential field and the potential field due to
			% the inclusion.% freq: the working frequency, a scalar
			%
			% Inputs:
			% freq: the working frequency, a scalar
			% s: index of the source
			% z0: center of the mesh
			% width: width of the mesh
			% N: number of points by side
			%
			% Outputs:
			% F: potential field u
			% F_bg: potential field U
			% Sx, Sy: coordinates in x and y-axis of the rectangular region
			% mask: binary mask where 0 indicates position where the potential may be numerically
			% undefined
			
			if obj.nbIncls > 1
				error('Multiple inclusions are not supported in the current version!');
			else
				D = obj.D{1};
				
				epsilon = width/(N-1)/5;
				
				[SX, SY, mask] = obj.D{1}.boundary_off_mask(z0, width, N, epsilon);
				
				SXi = SX(find(mask)); SYi = SY(find(mask));
				Zi = [SXi(:) SYi(:)]';
				SXo = SX(find(1-mask)); SYo = SY(find(1-mask));
				Zo = [SXo(:) SYo(:)]';
				
				% Resolution of the system 5.11
				xs = obj.cfg.src(s);
				sol = PDE.Helmholtz_R2.solve_forward(D, freq, xs, obj.pmeb, obj.pmtt, obj.pmeb_bg, obj.pmtt_bg);
				vphi = sol(1:size(sol,1)/2, :);
				vpsi = sol(size(sol,1)/2+1:end, :);
				
				k0 = sqrt(obj.pmtt_bg*obj.pmeb_bg)*freq;
				umU = ops.SingleLayer_H.eval(k0, D, vpsi, Zo);
				U = tools.Helmholtz.plane_wave(k0, Zo, xs);
				uo = umU(:)+U(:);
				
				ks = sqrt(obj.pmtt*obj.pmeb)*freq;
				ui =  ops.SingleLayer_H.eval(ks, D, vphi, Zi);
				
				F=zeros(size(mask));
				F(find(mask)) = ui;
				F(find(1-mask)) = uo;
				F_bg = tools.Helmholtz.plane_wave(k0, [SX(:) SY(:)]', xs);
				F_bg = reshape(F_bg, size(mask));
			end
		end
		
		out = data_simulation(obj, freq)
		
		plot_field(obj, F, F_bg, SX, SY, nbLine, varargin)
		
		out = reconstruct_SCT(obj, MSR, freq, ord, maxiter, tol)
		out = reconstruct_SCT_analytic(obj, MSR, freq, ord)
		
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
		
		function M = system_matrix_block_A11(freq, D, type, step, pmeb, pmtt)
			k = tools.Helmholtz.wavenb(freq, pmeb, pmtt);
			SH = ops.SingleLayer_H(k, D, type, step);
			M = SH.stiffmat;
		end
		
		function M = system_matrix_block_A12(freq, D, type, step, pmeb_bg, pmtt_bg)
			k0 = tools.Helmholtz.wavenb(freq, pmeb_bg, pmtt_bg);
			SH = ops.SingleLayer_H(k0, D, type, step);
			M = -1 * SH.stiffmat;
		end
		
		function M = system_matrix_block_A21(freq, D, type, step, pmeb, pmtt)
			k = tools.Helmholtz.wavenb(freq, pmeb, pmtt);
			Id = ops.Ident(D, type, step);
			KsH = ops.Kstar_H(k, D, type, step);
			M = (-1/2*Id.stiffmat + KsH.stiffmat)/pmeb;
		end
		
		function M = system_matrix_block_A22(freq, D, type, step, pmeb_bg, pmtt_bg)
			k0 = tools.Helmholtz.wavenb(freq, pmeb_bg, pmtt_bg);
			Id = ops.Ident(D, type, step);
			KsH = ops.Kstar_H(k0, D, type, step);
			M = - (1/2*Id.stiffmat + KsH.stiffmat)/pmeb_bg;
		end
		
		function [b] = source_vector(D, freq, Xi, pmeb_bg, pmtt_bg)
			% Compute the vector of right hand side in eq. 5.11 of [1], namely, U and \partialU\partial\nu,
			% where U is the plane wave: e^{k*i*<\xi,x>}
			
			k0 = tools.Helmholtz.wavenb(freq, pmeb_bg, pmtt_bg);
			s1 = tools.Helmholtz.plane_wave(k0, D.points, Xi);
			s2 = tools.Helmholtz.plane_wave_Dn(k0, D, Xi)/pmeb_bg;
			
			b = [s1; s2]; % Concatenuation gives a 2N-by-M, with N=D.nbPoints, M=number of
			% source directions
		end
		
		function [b] = source_vector_cw(D, freq, m, pmeb_bg, pmtt_bg)
			% Compute the vector of right hand side in eq. 5.11 of [1], namely, U and \partialU\partial\nu,
			% where U is the cylindrical wave: Jm(k_0|x|)e^{im theta_x}
			
			k0 = tools.Helmholtz.wavenb(freq, pmeb_bg, pmtt_bg);
			s1 = tools.Helmholtz.cylind_wave(k0, m, D.points);
			s2 = tools.Helmholtz.cylind_wave_Dn(k0, m, D)/pmeb_bg;
			
			b = [s1; s2];
		end
		
		function [sol] = solve_forward(D, freq, Xi, pmeb, pmtt, pmeb_bg, pmtt_bg)
			% Solve the forward problem 5.11 of [1] using plane wave as sources U.
			
			if length(freq)>1
				error('Frequency must be a scalar.');
			end
			
			% Construct the blocks of the system matrix
			matrix_A = PDE.Helmholtz_R2.system_matrix_block_A11(freq, D, 'P0', 1, pmeb, pmtt);
			matrix_B = PDE.Helmholtz_R2.system_matrix_block_A12(freq, D, 'P0', 1, pmeb_bg, pmtt_bg);
			matrix_C = PDE.Helmholtz_R2.system_matrix_block_A21(freq, D, 'P0', 1, pmeb, pmtt);
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
		
		function [sol] = solve_forward_cw(D, freq, m, pmeb, pmtt, pmeb_bg, pmtt_bg)
			% Solve the forward problem 5.11 of [1] using cylindrical wave
			% Um(x)=Jm(k_0|x|)e^{imtheta_x} as source.
			
			% Construct the blocks of the system matrix
			matrix_A = PDE.Helmholtz_R2.system_matrix_block_A11(freq, D, 'P0', 1, pmeb, pmtt);
			matrix_B = PDE.Helmholtz_R2.system_matrix_block_A12(freq, D, 'P0', 1, pmeb_bg, pmtt_bg);
			matrix_C = PDE.Helmholtz_R2.system_matrix_block_A21(freq, D, 'P0', 1, pmeb, pmtt);
			matrix_D = PDE.Helmholtz_R2.system_matrix_block_A22(freq, D, 'P0', 1, pmeb_bg, pmtt_bg);
			
			b = PDE.Helmholtz_R2.source_vector_cw(D, freq, m, pmeb_bg, pmtt_bg);
			
			matrix_BEM = [[matrix_A matrix_B] ; [matrix_C matrix_D]] ;
			sol = matrix_BEM\b;
		end
		
		function [F, F_bg] = evaluate_field(D, freq, xs, pmeb, pmtt, pmeb_bg, pmtt_bg, SX, SY, mask)
			% Evaluate the solution field u and U for a given source xs at the points of coordinates given by
			% SX and SY, which should not contain the interior of D.
			
			SXi = SX(find(mask)); SYi = SY(find(mask));
			Zi = [SXi(:) SYi(:)]';
			SXo = SX(find(1-mask)); SYo = SY(find(1-mask));
			Zo = [SXo(:) SYo(:)]';
			
			% Resolution of the system 5.11
			sol = PDE.Helmholtz_R2.solve_forward(D, freq, xs, pmeb, pmtt, pmeb_bg, pmtt_bg);
			vphi = sol(1:size(sol,1)/2, :);
			vpsi = sol(size(sol,1)/2+1:end, :);
			
			k0 = sqrt(pmtt_bg*pmeb_bg)*freq;
			umU = ops.SingleLayer_H.eval(k0, D, vpsi, Zo);
			U = tools.Helmholtz.plane_wave(k0, Zo, xs);
			uo = umU(:)+U(:);
			
			ks = sqrt(pmtt*pmeb)*freq;
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
