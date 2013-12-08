function out = data_simulation(obj, freq)
% Generate the multi-static response (MSR) matrix of the Helmholtz equation (the perturbed field
% u-U) by solving the linear system (5.11) of [1]

    if obj.nbIncls > 1
        error('NotImplemented: data_simulation for multi-inclusions.')
    else   
        if nargin < 2
            freq = obj.freq;
        end
        
        % rcv = obj.cfg.all_rcv;
        % src = obj.cfg.all_src;

        D = obj.D{1};
        
        for n=1:length(freq)
            MSR = zeros(obj.cfg.Nr, obj.cfg.Ns);
            
            % Resolution of the integral equation 5.11
            sol = PDE.Helmholtz_R2.solve_forward(D, freq(n), obj.cfg.all_src, obj.pmeb_bg, ...
                                                 obj.pmtt_bg);            
            vphi = sol(1:size(sol,1)/2, :);
            vpsi = sol(size(sol,1)/2+1:end, :);
            
            k_0 = sqrt(obj.pmtt_bg*obj.pmeb_bg)*freq(n);

            for s = 1:obj.cfg.Ns          
                MSR(:,s) = ops.SingleLayer_H.eval(k_0, D, vpsi(:,s), obj.cfg.rcv(s));
            end
            
            out.MSR{n} = transpose(MSR); % take the transpose, each row needs to correspond to a source
            out.vphi{n} = vphi;
            out.vpsi{n} = vpsi;
        end
        
        out.freq = freq;
    end
end
