classdef PulseImaging_R2 < PDE.Conductivity_R2
	% Class for the pulse imaging problem in the free space
	%
	% Reference:  Time-domain multiscale shape identification in electro-sensing using pulse-type sources.
	
	properties(Access=protected)
		% Relative to the linear system A[phi] = b.
		% A is a block matrix: (lambda I - K_{D_n}^*) on the diagonal, dSdn on the off diagonal
		% b can be the normal derivative of either the polynomials (for computation of GPTs), or
		% the Green's function (for simulation of data)
		% phi is the solution and it is used in either the computation of GPTs or the simulation
		% of MSR: u-G = \sum_l dSdn[phi_l]
		%
		% In the matrix A only the diagonal terms depend on the time step dt, the
		% other time-independent coefficients are stored in KsdS.
		% We reuse the variables KsdS and dGdn defined in the class |Conductivity_R2|
		
		%         KsdS=[] % system block matrix
		%         dGdn=[] % the normal derivative of the Green's function
		
		% dSLdn = {} % precomputed normal derivative of single layer operators
		waveform0 % pulse waveform
	end
	
	properties(SetAccess=protected)
		Ntime % length of waveform
		time % time interval on which the waveform is defined, time = [0,1,...Ntime]*dt
		dt % time step
	end
	
	
	%%% Recharge of some methods already defined in |Conductivity_R2| %%%
	
	methods(Access=protected) % Auxiliary functions
		
		function val = compute_dGdn(obj, sidx)
			% !!!!
			%  This function is defined in Conductivity_R2 class but
			%  recharged here.!!!!
			%
			% Construct the vector dG/dn of the given source. The output is
			% a 3D matrix of dimension: nbPoints X nbIncls X length(sidx)
			%
			% This will be used later to construct the right hand vector (1+alpha*dt)*dG/dn.
			% If the source is not given, compute for all sources.
			
			if nargin<2
				sidx = 1:obj.cfg.Ns_total;
			end
			
			nbPoints = obj.D{1}.nbPoints;
			val = zeros(nbPoints, obj.nbIncls, length(sidx));
			
			for i=1:obj.nbIncls
				toto = zeros(length(sidx), nbPoints);
				for s=1:length(sidx)
					psrc = obj.cfg.neutSrc(sidx(s)); % the positions of diracs of this source satisfying the neutrality condition (see acq.Concerntric.m)
					toto(s, :) = reshape(obj.cfg.neutCoeff, 1, []) * tools.Laplacian.Green2D_Dn(psrc, obj.D{i}.points, obj.D{i}.normal);
				end
				val(:, i, :) = toto';
			end
			
			% Equivalent, but more basic
			%             val = zeros(nbPoints*obj.nbIncls, length(sidx));
			%             idx=0;
			%
			%             for i=1:obj.nbIncls
			%                 toto = zeros(length(sidx), nbPoints);
			%                 for s=1:length(sidx)
			%                     psrc = obj.cfg.neutSrc(sidx(s)); % the positions of diracs of this source satisfying the neutrality condition (see acq.Concerntric.m)
			%                     toto(s, :) = reshape(obj.cfg.neutCoeff, 1, []) * tools.Laplacian.Green2D_Dn(psrc, obj.D{i}.points, obj.D{i}.normal);
			%                 end
			%                 val(idx+1:idx+nbPoints, :) = toto';
			%                 idx = idx+nbPoints;
			%             end
			% val = reshape(val, nbPoints, obj.nbIncls, []);
		end
		
		function Phi = compute_phi(obj, Ntime, s)
			% Compute the boundary function phi given the time interval and the source index
			% The function phi_s (s is the index of source) is the solution to
			% A[phi_s] = b_s
			% A is independent of s, we solve for multiple s by
			% [phi_1, phi_2, ...] = A\[b_1, b_2, ...]
			% Inputs:
			% Ntime: end time index of the time interval [0...Ntime]*dt
			% s: (optional) source indexes. If not given, compute for all
			% sources.
			%
			% Output:
			% Phi: Phi{t} is a 3D array of dimension (nbPoints X nbIncls X
			% nbSrc) at the time t. The meanings of these three dimensions are: 1. function values on the boundary, 2.inclusions, 3.sources
			
			nbPoints = obj.D{1}.nbPoints; % all inclusions have the same number of discretization points
			lambda = (obj.pmtt/obj.dt + obj.cnd + 1)./(obj.pmtt/obj.dt + obj.cnd - 1)/2;
			
			% Generate the system matrix
			Amat = asymp.CGPT.make_system_matrix_fast(obj.KsdS, lambda);
			[~, Amat_hlf]= asymp.CGPT.make_system_matrix_fast(obj.KsdS, 1/2*ones(1,obj.nbIncls));
			
			if nargin<3 % If the source index is not given, we solve for all sources
				dGdn = obj.dGdn;
			else
				dGdn = obj.dGdn(:,:,s);
			end
			
			nbSrc = size(dGdn,3);
			
			Phi{1} = zeros(nbPoints, obj.nbIncls, nbSrc);  % Initialization
			% Cm: the constant alpha_l / (alpha_l + dT)
			alphal = obj.pmtt ./ (obj.cnd-1);
			Cm = alphal ./ (alphal + obj.dt);
			
			for t=2:Ntime
				% Construct the right hand vector RHS
				W = zeros(nbPoints, obj.nbIncls, nbSrc);
				
				for i=1:obj.nbIncls
					V = zeros(nbPoints, nbSrc);
					for j=1:obj.nbIncls
						if j~=i
							% dSLdn = ops.dSLdn(obj.D{j}, 'P0', 1, obj.D{i}, 'P0', 1);
							V = V + Amat_hlf{i,j} * (squeeze(Phi{t-1}(:,j,:))); % for all possible sources, this is the same as apply -1*ops.dSLdn.fwd()
						end
					end
					
					W(:, i, :) = Cm(i) * (Amat_hlf{i,i} * (squeeze(Phi{t-1}(:,i,:))) + V - squeeze(dGdn(:, i, :)) * obj.waveform(t-1));
				end
				
				RHS = reshape(dGdn*obj.waveform(t) + W, nbPoints*obj.nbIncls, []);
				
				% function phi of all sources and inclusions
				Phi{t} = reshape(Amat\RHS, nbPoints, obj.nbIncls, []);
			end
		end
	end
	
	methods
		function obj = PulseImaging_R2(D, cnd, pmtt, waveform, dt, cfg)
			obj = obj@PDE.Conductivity_R2(D, cnd, pmtt, cfg);
			
			if ~isa(cfg, 'acq.Concentric')
				error('The configuration of the acquisition system must be object of the class acq.Concentric!');
			end
			
			if  cfg.nbDirac == 1
				warning('The source must fulfill the neutrality condition!');
			end
			
			if length(cnd)<obj.nbIncls || length(pmtt)<obj.nbIncls
				error('The value of conductivity and permittivity must be specified for each inclusion');
			end
			
			for n=1:obj.nbIncls
				if cnd(n)==1 || cnd(n)<0
					error('The conductivity constant must be positive and different from 1');
				end
			end
			
			% obj.cnd = cnd; obj.pmtt = pmtt;
			
			obj.waveform0 = reshape(waveform,1,[]);
			
			obj.dt = dt;
			
			% KsdS and dGdn have already been computed in Conductivity_R2. But dGdn needs to be re-computed.
			% This is not very clever but there is no better solution for the moment.
			% obj.KsdS = asymp.CGPT.make_block_matrix(obj.D);
			obj.dGdn = obj.compute_dGdn();
			
			%             obj.dSLdn = cell(obj.nbIncls); % cell of operators
			%             for i=1:obj.nbIncls
			%                 for j=1:obj.nbIncls
			%                     if j~=i
			%                         obj.dSLdn{i,j} = ops.dSLdn(obj.D{j}, 'P0', 1, obj.D{i}, 'P0', 1);
			%                     end
			%                 end
			%             end
		end
		
		function val = waveform(obj, t)
			if t <= 1 || t >= obj.Ntime
				val = 0;
			else
				val = obj.waveform0(t);
			end
		end
		
		function val = get.time(obj)
			val = obj.dt * (0:obj.Ntime-1);
		end
		
		function val = get.Ntime(obj)
			val = length(obj.waveform0);
		end
		
		out = data_simulation(obj, Tmax)
		
		% Calculate and plot the potential fields
		[F, F_bg, Sx, Sy, mask] = calculate_field(obj, Ntime, s, z0, width, N)
		% plot_field(obj, s, F, F_bg, Sx, Sy, nbLine, varargin)
		
		% Add white noise and reconstruct time-dependent CGPT
		CGPTt = addnoise_recon(obj, data, nlvl, ord, maxiter, tol, symmode, method, op)
	end
	
	methods(Static) % Utility functions
		
		function out = add_white_noise_global(data, nlvl)
			% Add white noise to simulated data. The noise is fixed by its
			% level (proportional to the data) and is treated globally for
			% all time instants in data, unlike add_white_noise() which
			% treats each instant independantly.
			% Input:
			% data: same as in add_white_noise()
			% nlvl: noise level
			
			out = data;
			Nt = length(data.MSR);
			
			toto = 0;
			
			for t=1:Nt
				toto = toto + norm(data.MSR{t}, 'fro')^2;
			end
			
			for t=1:Nt
				noise = randn(size(data.MSR{1}));
				out.MSR_noisy{t} = data.MSR{t} + noise * sqrt(toto/Nt/numel(data.MSR{1})) * nlvl;
			end
			
			data.sigma = sqrt(toto/Nt/numel(data.MSR{1})) * nlvl;
		end
		
	end
end
