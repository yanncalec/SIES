function out = data_simulation(obj, freq)
% Simulation of the data of electric fish
%
% Input:
% freq: a list of working frequencies
%
% Output:
% The output is a structure containing multiple terms. The variables whose name starting by
% 'v' are coefficient vectors of boundary element basis, and those starting by 'f' are function
% values. Those starting by a capital letter are measurements.
%
% We use n as the index of the frequency, i the index of the inclusion.
%
% vpsi{n}: each column is a solution vector psi of the linear system (A.5) with a
% given source
% vphi{n,i}: same as vpsi but for the vector phi
% fpsi{n}, fphi{n}: their function values.
%       fpsi{n} is 2D matrix: (nbBEM1, Ns_total), and fphi{n} is a 3D matrix: (nbBEM2, Ns_total, nbIncls)
% vpsi_bg, fpsi_bg: solution vectors and function values for A.2
% fpp{n}: the function (4.7): (1/2 - K^* - xi dDdn)(dudn - dUdn)
%
% Current{n}: surface current dudn
% Current_bg: surface current dUdn
%
% MSR{n}: multi-static response matrix for CGPT reconstruction, defined in the reference [2]
% SFR{n}: Space-Frequency response matrix dudn - dUdn, corresponded to the dipolar
% expansion, defined in section 4.2 of the reference [1]
% PP_SFR{n}: Post-processed response matrix, defined by eq (4.7). The negative sign before xi in
% eq (4.7) of [1] was a typo: it should be +.
%
% Convention (although a little bit strange) :
% For the matrix 'fxxx', each column corresponds to a source.
% For the measurement matrices (name starts by a capital letter), each row corresponds to a source.


%% Solve the forward problem for all positions

% SF,SU,PP below are coefficient vectors under boundary element basis
%
% 3D solution matrix of Eq (A.5), with inclusion, frequency dependent;
SF = zeros(obj.nbBEM1 + obj.nbIncls*obj.nbBEM2, length(freq), obj.cfg.Ns_total);
% Solution matrix of Eq (A.2), without inclusion
SU = zeros(obj.nbBEM1, obj.cfg.Ns_total);
% Post-processing
PP = zeros(obj.nbBEM1, length(freq), obj.cfg.Ns_total);

% Psi_int and Phi_int will be used later to reinforce the zero mean constraint of the solution
Psi_int = obj.Psi'*obj.Omega.sigma(:); % boundary integral of P1 elements \int_Omega psi
Phi_int = zeros(obj.nbBEM2, obj.nbIncls) ;
for i=1:obj.nbIncls
	Phi_int(:,i) = obj.D{i}.sigma(:); % boundary integral of P0 elements \int_Omega phi
end

% The first for loop is on the fish's position (different Omega), because it is more expensive to
% build the block matrices depending on Omega (A, B, C). Remark that both the system matrix and the
% right hand vector are Omega-dependent.
for s=1:obj.cfg.Ns_total
	% Construct the blocks of the system matrix
	Omega = obj.cfg.Bodies(s); % Fish's body at s-th position
	matrix_A = PDE.ElectricFish.system_matrix_block_A(Omega, obj.typeBEM1, obj.stepBEM1, obj.impd);
	matrix_B = PDE.ElectricFish.system_matrix_block_B(Omega, obj.typeBEM1, obj.stepBEM1, obj.D, obj.typeBEM2, obj.stepBEM2);
	matrix_C = PDE.ElectricFish.system_matrix_block_C(Omega, obj.typeBEM1, obj.stepBEM1, obj.D, obj.typeBEM2, obj.stepBEM2, obj.impd);
	
	% Resolution of system A.2 - background field equation
	matrix_BEM = [matrix_A; reshape(Psi_int,1,[])];
	SU(:,s) = matrix_BEM\[obj.dHdn(1:obj.nbBEM1,s); 0];
	
	% The second for loop is on the frequency. We build then the frequency-dependent system
	% matrix and solve it
	for n=1:length(freq)
		lambda = asymp.CGPT.lambda(obj.cnd, obj.pmtt, freq(n));
		
		matrix_D = asymp.CGPT.make_system_matrix_fast(obj.KsdS, lambda);
		
		% Generate the system matrix and resolution of system A.5
		matrix_BEM = [[matrix_A matrix_B] ; [matrix_C matrix_D]; [reshape(Psi_int, 1, []) ...
			zeros(1,obj.nbBEM2*obj.nbIncls)]; [zeros(1,obj.nbBEM1) reshape(Phi_int,1,[])]];
		SF(:,n,s) = matrix_BEM\[obj.dHdn(:,s); 0; 0]; % function (psi,phi) of the source s and all inclusions
		
		% Post-processing (1/2-K^*+xi dDdn)[dudn - dUdn]
		%
		% We calculate not its function value but the representation
		% coefficients under P1-basis. For this, first we get the stiffness
		% matrix of the operator, then apply it on the synthesis coefficient of
		% [dudn - dUdn], which yields in turn the analysis coefficient.
		% Finally, we apply the inverse of Gram matrix on the analysis
		% coefficient {a_n} to get the synthesis coefficient {s_n} :
		%
		% f = \sum_n s_n h_n, and a_n = <f,h_n>, with h_n a basis, then s = G^-1 a
		
		% Using the fact that dudn|+=psi_u (psi_u is the function psi of A.5),
		% and dUdn|+=psi_U (psi_U is the function psi of A.2), we see that
		% Postprocessing(du|+ - dU|+) = Postprocessing(psi - psiU) = S_D[phi]
		%
		% The Gram matrix can not be removed otherwise the dipolar
		% expansion will not match the data
		PP(:,n,s) = -obj.Grammatrix\(matrix_B * SF(obj.nbBEM1+1:end, n, s)); % minus sign because matrix_B = -S_D
	end
end

out.freq = freq; % save the frequency list

%% Arranging the solutions
out.vpsi_bg = SU; % P1 basis coefficient of the solution psi to the system (A.2)
out.fpsi_bg = tools.BEM.interpolation(obj.Psi, SU, []); % its function value

% vphi and fphi are the same thing.
out.vphi=cell(length(freq), obj.nbIncls);

for n = 1:length(freq)
	sol = squeeze(SF(:, n, :)); % the space solution corresponding to the n-th frequency
	out.vpsi{n} = sol(1:obj.nbBEM1, :); % vector psi, see (eq A.5)
	% interpolation to convert the coefficient vector to function value
	out.fpsi{n} = tools.BEM.interpolation(obj.Psi, out.vpsi{n}, []);
	
	idx=obj.nbBEM1;
	out.fphi{n} = zeros(obj.nbBEM2, obj.cfg.Ns_total, obj.nbIncls);
	for i=1:obj.nbIncls
		out.vphi{n,i} = sol(idx+1:idx+obj.nbBEM2, :); % vector phi, see (eq A.5)
		idx = idx+obj.nbBEM2;
		out.fphi{n}(:,:,i) = out.vphi{n,i}; % no need of interpolation for P0 basis since the step = 1
	end
	
	% Post-proccessed SFR evaluated on the whole body
	vpp = squeeze(PP(:, n, :)); % Remind that PP is synthesis coefficient of the P1 basis
	out.fpp{n} =  tools.BEM.interpolation(obj.Psi, vpp, []);
	
	% Compute the MSR matrix
	out.MSR{n} = zeros(obj.cfg.Ns_total, obj.cfg.Nr);
	
	for i=1:obj.nbIncls
		toto = zeros(obj.cfg.Ns_total, obj.cfg.Nr);
		for s=1:obj.cfg.Ns_total
			% Evaluate the single layer potential S_D[phi] on fish's body
			rcv = obj.cfg.rcv(s); % receivers corresponding to the s-th source
			toto(s,:) = ops.SingleLayer.eval(obj.D{i}, out.vphi{n,i}(:,s), rcv);
		end
		out.MSR{n} = out.MSR{n}+toto; % active receivers have already been automatically used
	end
end

%% For ease of reference, rename some variables.

% The measurement is dUdn|+, which equals to psi.
Current_bg = transpose(out.fpsi_bg);

% Take the measurement at activate receptors.
for n = 1:length(freq)
	Current = transpose(out.fpsi{n}); % dudn
	SFR = Current - Current_bg; % dudn - dUdn
	PP_SFR = transpose(out.fpp{n}); % (1/2-K^*+xi dDdn)[dudn - dUdn]
	
	out.Current{n} = Current(:, obj.cfg.idxRcv);
	out.SFR{n} = SFR(:, obj.cfg.idxRcv);
	out.PP_SFR{n} = PP_SFR(:, obj.cfg.idxRcv);
end

out.Current_bg = Current_bg(:, obj.cfg.idxRcv);
