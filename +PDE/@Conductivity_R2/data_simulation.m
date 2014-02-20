function out = data_simulation(obj, freq)
% Simulation of the MSR matrix:
% MSR(s,r) = u-Hs(xr) = sum_l S_{D_l}[phi_l^s](x_r)
% Inputs:
% freq: a list of working frequencies
% Output:
% out.MSR: a list of MSR matrix of all working frequencies. The coefficient at the s-th
% row and r-th column is the data obtained with the s-th source and the r-th that respond
% to this source.

    if nargin<2
        freq = 0;
    end

    for f=1:length(freq)
        Phi = obj.compute_phi(freq(f));

        % Compute the MSR matrix by evaluating the single layer potential
        MSR = zeros(obj.cfg.Ns_total, obj.cfg.Nr);

        for i=1:obj.nbIncls
            toto = zeros(obj.cfg.Ns_total, obj.cfg.Nr);
            for s=1:obj.cfg.Ns_total
                rcv = obj.cfg.rcv(s); % receivers corresponding to the s-th source
                toto(s,:) = ops.SingleLayer.eval(obj.D{i}, Phi{i}(:,s), rcv);
            end
            MSR = MSR+toto;
        end                
        out.MSR{f} = MSR;                
    end
    out.freq = freq; % save the frequency list

end        
