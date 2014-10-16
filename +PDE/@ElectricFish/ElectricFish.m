classdef ElectricFish < PDE.Small_Inclusions
	% Class for modeling the weakly electric fish
	% Reference:
	% [1] Modeling active electrolocation in weakly electric fish. Ammari, Boulier, Garnier, version 28, June, 2012
	% [2] Shape identification in electrolocation, Ammari, Boulier, Garnier and Wang, PNAS 2014
	
	properties(Access=protected)
		% Similar to Conductivity_R2 class, the system matrix is build by blocks as
		% [A         B        ]
		% [C   lambda I - KsdS]
		
		KsdS=[] % same as KsdS of Conductivity_R2
		dHdn=[] % the normal derivative of H, or the right hand vector in A[phi]=b
	end
	
	properties(SetAccess = protected)
		Omega % domain for the fish's body, an object of C2boundary class
		impd % impedance of the skin, a scalar
		
		cnd % conductivity constant corresponding to each inclusion, an array
		pmtt % permittivity constant corresponding to each inclusion, an array
		
		typeBEM1 % type of boundary elements for the fish, 'P0', 'P1' etc.
		nbBEM1 % number of (P1) boundary elements for the fish's body
		stepBEM1 % downsampling factor of the boundary elements for the fish
		
		typeBEM2 % type of boundary elements for inclusions (same for all)
		nbBEM2 % number of boundary elements for inclusions
		stepBEM2 % downsampling factor of the boundary elements for inclusions
		
		Psi % P1 basis for fish's body
		Grammatrix % Gram matrix of Psi: Psi' * Psi
	end
	
	methods(Access=protected)
		% Compute the normal derivative of the function H in eq A.5, the right hand vector
		val  = compute_dHdn(obj, sidx)
	end
	
	methods
		function obj = ElectricFish(D, cnd, pmtt, cfg, stepBEM)
			% obj = ElectricFish(Omega, D, cfg, stepBEM)
			% Constructor of ElectriFish class
			% Inputs:
			% D: exterior small inclusions
			% cnd, pmtt: conductivity and permittivity constants of D
			% cfg: configuration acq.Fish_circle
			% stepBEM: sampling step for the P1 BEM
			
			obj = obj@PDE.Small_Inclusions(D, cfg);
			
			if ~isa(cfg, 'acq.Fish_circle')
				error('Configuration must be an object of acq.Fish_circle.');
			end
			
			% make a copy of Omega and impd which will be frequently used
			obj.Omega = cfg.Omega0;
			obj.impd = cfg.impd;
			
			obj.cnd = cnd; obj.pmtt = pmtt;
			
			obj.typeBEM1 = 'P1'; obj.typeBEM2 = 'P0';
			
			if mod(obj.Omega.nbPoints, stepBEM)
				error('Sampling step for the fish is invalid');
			end
			
			obj.stepBEM1 = stepBEM; obj.stepBEM2 = 1;
			
			% P1 elements for fish body;
			obj.Psi = tools.BEM.P1_basis(obj.Omega.nbPoints, obj.stepBEM1);
			
			obj.KsdS = asymp.CGPT.make_block_matrix(obj.D);
			obj.dHdn = obj.compute_dHdn();
		end
		
		function val = get.Grammatrix(obj)
			val = obj.Psi' * diag(obj.Omega.sigma) * obj.Psi;
		end
		
		function val = get.nbBEM1(obj)
			val = floor(obj.Omega.nbPoints / obj.stepBEM1);
		end
		
		function val = get.nbBEM2(obj)
			val = floor(obj.D{1}.nbPoints / obj.stepBEM2);
		end
		
		function plot(obj, idx, varargin)
			for n=1:obj.nbIncls
				plot(obj.D{n}, varargin{:}); hold on;
			end
			% plot then the acquisition system
			if nargin<2
				idx = [];
			end
			obj.cfg.plot(idx, varargin{:});
		end
		
		[F, F_bg, SX, SY, mask] = calculate_field(obj, f, s, z0, width, N, fpsi_bg, fpsi, fphi)
		plot_field(obj, s, F, F_bg, SX, SY, subfig, varargin)
		
		out = data_simulation(obj, freq)
		
		out = reconstruct_CGPT(obj, MSR, Current, ord, maxiter, tol, symmode)
		out = reconstruct_PT(obj, SFR, Current_bg, maxiter, tol, symmode)
		
		SF_MUSIC(obj, s, PP_SFR, SX, SY, varargin) % TODO
		
	end
	
	methods(Static)
		function matrix_A = system_matrix_block_A(Omega, type1, step1, impd)
			% Post-processing operator (1/2-K_Omega^*+xi dD_Omega dn)
			Id = ops.Ident(Omega, type1, step1);
			Ks = ops.Kstar(Omega, type1, step1);
			dDLdn = ops.dDLdn(Omega, type1, step1);
			matrix_A = Id.stiffmat/2 - Ks.stiffmat + impd * dDLdn.stiffmat;
		end
		
		function matrix_B = system_matrix_block_B(Omega, type1, step1, D, type2, step2)
			% Operator -dS_{D_l}dn restricted on Omega
			nbIncls = length(D);
			toto = cell(1, nbIncls);
			for n=1:nbIncls
				dSLdn = ops.dSLdn(D{n}, type2, step2, Omega, type1, step1);
				toto{n} = -1*dSLdn.stiffmat;
			end
			matrix_B = cell2mat(toto);
		end
		
		function matrix_C = system_matrix_block_C(Omega, type1, step1, D, type2, step2, impd)
			% Operator -dS_{Omega}dn + xi*dD_{Omega}dn restricted on D_l
			nbIncls = length(D);
			toto = cell(nbIncls,1);
			for n=1:nbIncls
				dSLdn = ops.dSLdn(Omega, type1, step1, D{n}, type2, step2);
				dDLdn = ops.dDLdn(Omega, type1, step1, D{n}, type2, step2);
				toto{n} = -1*dSLdn.stiffmat + impd*dDLdn.stiffmat;
			end
			matrix_C = cell2mat(toto);
		end
		
		function out = add_white_noise(data, nlvl)
			% Add white noise to simulated data, see the function data_simulation
			% nlvl: noise level
			
			out = data;
			mode = 1; % treat each source independently
			rowmajor = 1; % procede row-by-row
			
			% add real white noise
			for n=1:length(data.SFR)
				out.MSR_noisy{n} = tools.add_white_noise(real(out.MSR{n}), nlvl, mode, rowmajor);
				out.SFR_noisy{n} = tools.add_white_noise(real(out.SFR{n}), nlvl, mode, rowmajor);
				out.PP_SFR_noisy{n} = tools.add_white_noise(real(out.PP_SFR{n}), nlvl, mode, rowmajor);
				out.Current_noisy{n} = tools.add_white_noise(real(out.Current{n}), nlvl, mode, rowmajor);
			end
			
			% add imaginary white noise
			for n=1:length(data.SFR)
				out.MSR_noisy{n} = out.MSR_noisy{n} + 1i * tools.add_white_noise(imag(out.MSR{n}), nlvl, mode, rowmajor);
				out.SFR_noisy{n} = out.SFR_noisy{n} + 1i * tools.add_white_noise(imag(out.SFR{n}), nlvl, mode, rowmajor);
				out.PP_SFR_noisy{n} = out.PP_SFR_noisy{n} + 1i * tools.add_white_noise(imag(out.PP_SFR{n}), nlvl, mode, rowmajor);
				out.Current_noisy{n} = out.Current_noisy{n} + 1i * tools.add_white_noise(imag(out.Current{n}), nlvl, mode, rowmajor);
			end
			
			% only real noise for dUdn since it is real
			out.Current_bg_noisy = tools.add_white_noise(out.Current_bg, nlvl, mode, rowmajor);
		end
		
		% Linear operators involved in the reconstruction of CGPT and PT
		out = make_linop_CGPT(cfg, Current, impd, ord, symmode)
		out = make_linop_PT(cfg, Current_bg, impd, symmode)
	end
end

