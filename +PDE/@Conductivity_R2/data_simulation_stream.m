function out = data_simulation_stream(obj, Z, A, freq)
% Simulation of the MSR matrices of a mobile target.
%
% Inputs:
% Z: positions (trajectory) of the mobile target, dimension 2 X Ntime
% A: angular positions (rotation, orientation) of the mobile target, dimension 1 X Ntime
% freq: working frequency, 0 by default
%
% Output:
% out.MSR: a list of MSR matrix corresponding to the trajectory.

    if nargin<4
        freq = 0;
    else
        if iscell(freq)
            error('Data simulation of mobile target with Multi-frequency is not supported!');
        end
    end
    
    if obj.nbIncls>1
        error('Data simulation of multi-inclusions is not supported!');
    end
        
    if length(Z) ~= length(A)
        error('Position and orientation of the mobile target must have the same length!')
    else        
        ntime = length(Z);
    end

    Phi = obj.compute_phi(freq);
    MSR = {};
    
    for t=1:ntime
        Dt = (obj.D{1}<A(t)) + Z(:,t); % target after the rigid motion

        % Compute the MSR matrix by evaluating the single layer potential
        toto = zeros(obj.cfg.Ns_total, obj.cfg.Nr);

        for s=1:obj.cfg.Ns_total
            rcv = obj.cfg.rcv(s); % receivers corresponding to the s-th source
            toto(s,:) = ops.SingleLayer.eval(Dt, Phi{1}(:,s), rcv);
        end
        MSR{t} = toto;
    end
    
    % Convert the data from cell to matrix format
    MSR_mat = zeros(numel(toto), ntime);
    for t=1:ntime
        MSR_mat(:,t) = reshape(MSR{t}, [], 1);
    end
    
    out.MSR = MSR;
    out.MSR_mat = MSR_mat;
    
    out.freq = freq;
    out.ntime = ntime;
    
    out.Z = Z;
    out.A = A;
    
    
end
