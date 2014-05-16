function out = data_simulation(obj, Ntime)
% Simulation of the time-dependent MSR matrix:
% MSR(s,r) = (u-U)(xr) = sum_l S_{D_l}[phi_l^s](x_r)
%
% Inputs:
% Ntime: end time index (if not given then simulate for the whole time
% duration of the pulse signal)
%
% Output:
% out.MSR: a list of MSR matrix stream. The coefficient at the s-th
% row and r-th column is the data obtained with the s-th source and the r-th that respond
% to this source.

    if nargin < 2
        Ntime = obj.Ntime;
    end
    
    Phi = obj.compute_phi(Ntime); % Compute the function phi for all sources 
    out.MSR = {};
    out.Phi = Phi;
    
    for t=1:Ntime
        % Compute the MSR matrix by evaluating the single layer potential
        MSR = zeros(obj.cfg.Ns_total, obj.cfg.Nr);

        for i=1:obj.nbIncls
            toto = zeros(obj.cfg.Ns_total, obj.cfg.Nr);
            for s=1:obj.cfg.Ns_total
                rcv = obj.cfg.rcv(s); % receivers corresponding to the s-th source

                toto(s,:) = ops.SingleLayer.eval(obj.D{i}, Phi{t}(:,i,s), rcv);
            end
            MSR = MSR+toto;
        end                
        out.MSR{t} = MSR;                
    end
    out.Ntime = Ntime; % save the time list    
end        
