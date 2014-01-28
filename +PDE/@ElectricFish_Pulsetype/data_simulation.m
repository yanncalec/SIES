function out = data_simulation(obj, Nt)
% Simulation of the data of electric fish
%
% Input:
% Nt: number of time frames
%
% Output:
% The output is a structure containing multiple terms. The variables whose name starting by
% 'v' are coefficient vectors of boundary element basis, and those starting by 'f' are function
% values. Those starting by a capital letter are measurements.
%
% We use n as the index of the time, i the index of the inclusion.
%
% vpsi{n}: each column is a solution vector psi of the linear system with a given source
% vphi{n,i}: same as vpsi but for the vector phi
% fpsi{n}, fphi{n}: their function values.
%       fpsi{n} is 2D matrix: (nbBEM1, Ns_total), and fphi{n} is a 3D matrix: (nbBEM2, Ns_total, nbIncls)
% vpsi_bg, fpsi_bg: solution vectors and function values for 
% fpp{n}: the function: (1/2 - K^* + xi dDdn)(dudn - dUdn)
%
% Current{n}: surface current dudn
% Current_bg: surface current dUdn
%
% MSR{n}: multi-static response matrix for CGPT reconstruction, defined in the reference [2]
% STR{n}: Space-Time response matrix dudn - dUdn
% PP_STR{n}: Post-processed response matrix.
%
% Convention (although a little bit strange) :
% For the matrix 'fxxx', each column corresponds to a source.
% For the measurement matrices (name starts by a capital letter), each row corresponds to a source.


%% Solve the forward problem for all positions and all times
%
% ST,SU,PP below are coefficient vectors under boundary element basis

% 3D solution matrix of psi and phi, with inclusion and time dependent
ST = zeros(obj.nbBEM1 + obj.nbIncls*obj.nbBEM2, Nt, obj.cfg.Ns_total);
% Solution matrix, without inclusion (separable so there is no time here)
SU = zeros(obj.nbBEM1, obj.cfg.Ns_total);
% Post-processing
PP = zeros(obj.nbBEM1, Nt, obj.cfg.Ns_total);

% Psi_int and Phi_int will be used later to reinforce the zero mean constraint of the solution
Psi_int = obj.Psi'*obj.Omega.sigma(:); % boundary integral of P1 elements \int_Omega psi
Phi_int = zeros(obj.nbBEM2, obj.nbIncls) ;
for i=1:obj.nbIncls
    Phi_int(:,i) = obj.D{i}.sigma(:); % boundary integral of P0 elements \int_Omega phi
end

% The first for loop is on the fish's position (different Omega), because it is more expensive to
% build the block matrices depending on Omega (A, B, C). Remark that both the system matrix (block A,B,C) and the
% right hand vector are Omega-dependent.

% Block matrix D is independent of Omega
alpha = obj.pmtt./(obj.cnd-1);
lambda = (obj.cnd + 1)./(obj.cnd - 1)/2 + obj.time_step./((obj.time_step + alpha).*(obj.cnd - 1));
KsdS = asymp.CGPT.make_block_matrix(obj.D);
matrix_D = asymp.CGPT.make_system_matrix_fast(KsdS, lambda);
matrix_Dtilde = asymp.CGPT.make_system_matrix_fast(KsdS, 1/2*ones(1,obj.nbIncls));

for s=1:obj.cfg.Ns_total
    % Construct the blocks of the system matrix
    Omega = obj.cfg.Bodies(s); % Fish's body at s-th position
    
    matrix_A = PDE.ElectricFish.system_matrix_block_A(Omega, obj.typeBEM1, obj.stepBEM1, obj.impd);
    matrix_B = PDE.ElectricFish.system_matrix_block_B(Omega, obj.typeBEM1, obj.stepBEM1, obj.D, obj.typeBEM2, obj.stepBEM2);
    matrix_C = PDE.ElectricFish.system_matrix_block_C(Omega, obj.typeBEM1, obj.stepBEM1, obj.D, obj.typeBEM2, obj.stepBEM2, obj.impd);
    
    % Resolution of the background field equation
    matrix_BEM_U = [matrix_A; reshape(Psi_int,1,[])];
    SU(:, s) = matrix_BEM_U\[obj.dHdn(1:obj.nbBEM1,s); 0];
    
    % The system matrix is invariant to time and depends only on the source
    matrix_BEM = [[matrix_A matrix_B] ; [matrix_C matrix_D]; [reshape(Psi_int, 1, []) ...
        zeros(1,obj.nbBEM2*obj.nbIncls)]; [zeros(1,obj.nbBEM1) reshape(Phi_int,1,[])]];
    
    % We suppose that at time 0, the wave profile h(0)=0. This
    % implies that the solution vector is also 0.    
    bvec = obj.dHdn(:, s)*obj.waveform(1);

    % The second for loop is on the time.    
    for n=2:Nt        
        % psi and phi of time n-1
        psit = squeeze(ST(1:obj.nbBEM1, n-1, s));
        phit = squeeze(ST((obj.nbBEM1+1):end, n-1, s));
        
        % innovation
        bvec_D = bvec(obj.nbBEM1+1:end);
        toto = matrix_C * psit + matrix_Dtilde * phit - bvec_D;
        innov = kron(diag(alpha./(obj.time_step+alpha)), eye(obj.nbBEM2)) * toto;
        
        % b vector of time n
        bvec = obj.dHdn(:, s)*obj.waveform(n); % dH/dn at the time n
        rhs_vec = bvec + [zeros(obj.nbBEM1, 1); innov]; % the right hand side vector at time n
        
        ST(:, n, s) = matrix_BEM\[rhs_vec; 0; 0]; % function (psi,phi) of the source s and all inclusions
        
        % Post-processing (1/2-K^*+xi dDdn)[dudn - dUdn], same as in
        % ElectricFish case
        PP(:, n, s) = -obj.Grammatrix\(matrix_B * ST(obj.nbBEM1+1:end, n, s));
    end
end

%% Arranging the solutions
out.vpsi_bg = SU; % P1 basis coefficient of the solution psi to the system (A.2)
out.fpsi_bg = tools.BEM.interpolation(obj.Psi, SU, []); % its function value

% vphi and fphi are the same thing.
out.vphi=cell(Nt, obj.nbIncls);

for n = 1:Nt
    sol = squeeze(ST(:, n, :)); % the space solution corresponding to the n-th time frame
    out.vpsi{n} = sol(1:obj.nbBEM1, :); % vector psi
    % interpolation to convert the coefficient vector to function value
    out.fpsi{n} = tools.BEM.interpolation(obj.Psi, out.vpsi{n}, []);
    
    idx=obj.nbBEM1;
    out.fphi{n} = zeros(obj.nbBEM2, obj.cfg.Ns_total, obj.nbIncls);
    for i=1:obj.nbIncls
        out.vphi{n,i} = sol(idx+1:idx+obj.nbBEM2, :); % vector phi
        idx = idx+obj.nbBEM2;
        out.fphi{n}(:,:,i) = out.vphi{n,i}; % no need of interpolation for P0 basis since the step = 1
    end
    
    % Post-proccessed STR evaluated on the whole body
    vpp = squeeze(PP(:, n, :)); % Remind that PP is synthesis coefficient of the P1 basis
    out.fpp{n} =  tools.BEM.interpolation(obj.Psi, vpp, []);
    
    % Compute the MSR matrix
    out.MSR{n} = zeros(obj.cfg.Ns_total, obj.cfg.Nr);
    
    for i=1:obj.nbIncls
        toto = zeros(obj.cfg.Ns_total, obj.cfg.Nr);
        for s=1:obj.cfg.Ns_total
            % Evaluate the single layer potential S_D[phi] on fish's body, recall that phi =
            % (lambda I - K_D^*)^-1 [dH/dn]
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
for n = 1:Nt
    Current = transpose(out.fpsi{n}); % dudn
    SFR = Current - Current_bg; % dudn - dUdn
    PP_SFR = transpose(out.fpp{n}); % (1/2-K^*+xi dDdn)[dudn - dUdn]
    
    out.Current{n} = Current(:, obj.cfg.idxRcv);
    out.SFR{n} = SFR(:, obj.cfg.idxRcv);
    out.PP_SFR{n} = PP_SFR(:, obj.cfg.idxRcv);
end

out.Current_bg = Current_bg(:, obj.cfg.idxRcv);
