classdef PulseImaging_R2 < PDE.Small_Inclusions
% Class for the pulse imaging problem in the free space

    properties(Access=protected)
        % Relative to the linear system A[phi] = b.
        % A is a block matrix: (lambda I - K_{D_n}^*) on the diagonal, dSdn on the off diagonal
        % b can be the normal derivative of either the polynomials (for computation of GPTs), or
        % the Green's function (for simulation of data)
        % phi is the solution and it is used in either the computation of GPTs or the simulation
        % of MSR: u-G = \sum_l dSdn[phi_l]
        %
        % In the matrix A only the diagonal terms depend on the time step dt (via lambda I), the
        % other time-independent coefficients are stored in KsdS.
        KsdS=[] % system block matrix
        dGdn=[] % the normal derivative of the Green's function

        waveform0 % pulse waveform
    end
    
    properties(SetAccess=protected)
        cnd % conductivity constant corresponding to each inclusion, an array
        pmtt % permittivity constant corresponding to each inclusion, an array

        Ntime % length of waveform
        time % time interval on which the waveform is defined, time = [0,1,...Ntime]*dt
        dt % time step        
    end
    
    methods(Access=protected) % Auxiliary functions
        function val = compute_dGdn(obj, sidx)
        % Construct the right hand vector of the given source.
        % If the source is not given, compute for all sources
            
            if nargin<2
                sidx = 1:obj.cfg.Ns_total;
            end
            
            nbPoints = obj.D{1}.nbPoints;
            val = zeros(nbPoints*obj.nbIncls, length(sidx)); 
            idx=0;

            for i=1:obj.nbIncls
                toto = zeros(length(sidx), nbPoints);
                for s=1:length(sidx)
                    psrc = obj.cfg.neutSrc(sidx(s));
                    toto(s, :) = reshape(obj.cfg.neutCoeff, 1, []) * tools.Laplacian.Green2D_Dn(psrc, obj.D{i}.points, obj.D{i}.normal);
                end
                val(idx+1:idx+nbPoints, :) = toto';
                idx = idx+nbPoints; 
            end
            %             % Or by
            %             val = zeros(nbPoints, obj.nbIncls, length(sidx));
            %             val(:, i, :) = toto';
            %             val = reshape(val, nbPoints*obj.nbIncls, []);
        end
        
        function Phi = compute_phi(obj, Ntime, s)
        % Compute the boundary function phi given the time interval and the source index
        % The function phi_s (s is the index of source) is the solution to 
        % A[phi_s] = b_s
        % A is independent of s, we solve for multiple s by
        % [phi_1, phi_2, ...] = A\[b_1, b_2, ...]
        % Inputs:
        % Ntime: end time index of the time interval [0...Ntime]*dt
        % s: (optional) source index
        %
        % Output:
        % Phi: Phi{t} is a 3D array of dimension (nbPoints X nbIncls X
        % nbSrc) at the time t.
        
            nbPoints = obj.D{1}.nbPoints;
            lambda = (obj.pmtt/obj.dt + obj.cnd + 1)./(obj.pmtt/obj.dt + obj.cnd - 1)/2;

            % Generate the system matrix
            Amat = asymp.CGPT.make_system_matrix_fast(obj.KsdS, lambda);
            [~, Amat_hlf]= asymp.CGPT.make_system_matrix_fast(obj.KsdS, 1/2*ones(1,obj.nbIncls));
            
            if nargin<3
                dGdn = obj.dGdn;
            else
                dGdn = obj.compute_dGdn(s);
            end
            nbSrc = size(dGdn,2);

            Phi{1} = zeros(nbPoints, obj.nbIncls, nbSrc);            
            Cm = 1./(1 + obj.dt./(obj.pmtt ./ (obj.cnd-1)));

            for t=2:Ntime
                % Construct the right hand vector b
                
                W = zeros(nbPoints*obj.nbIncls, nbSrc);
                idx = 0;
                
                for i=1:obj.nbIncls
                    V = zeros(nbPoints, nbSrc);
                    
                    for j=1:obj.nbIncls
                        if j~=i
                            dSLdn = ops.dSLdn(obj.D{j}, 'P0', 1, obj.D{i}, 'P0', 1);
                            V = V + dSLdn.fwd(squeeze(Phi{t-1}(:,j,:))); % for all possible sources, see the function Operators.fwd
                        end
                    end
                    
                    toto = Amat_hlf{i,i} * squeeze(Phi{t-1}(:,i,:)) - V - dGdn(idx+1:idx+nbPoints, :) * obj.waveform(t-1);
                    W(idx+1:idx+nbPoints, :) = Cm(i) * toto;
                end
                
                RHS = dGdn*obj.waveform(t) + W;
                
                phi = Amat\RHS; % function phi of all sources and inclusions                
                Phi{t} = reshape(phi, nbPoints, obj.nbIncls, []); % the three dimensions are 1. function values on the boundary, 2.inclusions, 3.sources                
            end
        end
    end
    
    methods
        function obj = PulseImaging_R2(D, cnd, pmtt, waveform, dt, cfg)
            obj = obj@PDE.Small_Inclusions(D, cfg);

            if ~isa(cfg, 'acq.Concentric')
                error('Configuration of acquisition system must be object of the class acq.Concentric!');
            end
            
            if  cfg.nbDirac == 1
                error('Source must fulfill the neutrality condition!');
            end
            
            if length(cnd)<obj.nbIncls || length(pmtt)<obj.nbIncls
                error('The value of conductivity and permittivity must be specified for each inclusion');
            end
            
            for n=1:obj.nbIncls
                if cnd(n)==1 || cnd(n)<0
                    error('The conductivity constant must be positive and different from 1');
                end
            end
            
            obj.cnd = cnd; obj.pmtt = pmtt;
            
            obj.waveform0 = reshape(waveform,1,[]);
            
            obj.dt = dt;
            
            obj.KsdS = asymp.CGPT.make_block_matrix(obj.D); 
            obj.dGdn = obj.compute_dGdn();
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
        [F, F_bg, Sx, Sy, mask] = calculate_field(obj, freq, s, z0, width, N)        
        plot_field(obj, s, F, F_bg, Sx, Sy, nbLine, varargin)
        
        % Reconstruction of contracted GPT matrix
        out = reconstruct_CGPT(obj, MSR, ord, maxiter, tol, symmode) 
        out = reconstruct_CGPT_analytic(obj, MSR, ord)
        
    end
        
    methods(Static) % Utility functions                            
        A = make_matrix_A(Xs, Z, order)
        out = make_linop_CGPT(cfg, ord, symmode)  % construct the linear operator - CGPT

        function [waveform, dt] = make_pulse(Tmax, Ntime)
        % Make pulse waveform
          
            x= sym('x');
            f=exp(-x^2*2); g=simplify(diff(f,5));
            h = inline(g);
            waveform = h(linspace(-Tmax,Tmax,Ntime));
            dt = 2*Tmax/Ntime;

            % Method 2: first order derivative of a Gaussian
            % Tmax = 0.1;
            % x = linspace(-Tmax/2,Tmax/2,Ntime);
            % sgm = Tmax^2;
            % y = [0 exp(-x.^2/(2*sgm^2))/(sqrt(2*pi)*sgm)];
            % dt = Tmax/Ntime;
            % waveform = diff(y)/dt;

            %             Tmax = 0.1; %[0, Tmax]
            %             hfunc_inline = inline('exp(-(x-m1).^2/2/s1^2)-exp(-(x-m2).^2/2/s2^2)');
            %             hfunc = @(Tmax,x)hfunc_inline(Tmax/2-Tmax*0.01,Tmax/2+Tmax*0.01,Tmax*0.1,Tmax*0.1, x);
            %             waveform = hfunc(Tmax, linspace(0,Tmax,99)); waveform = [0, waveform]*10;
        end
        
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
