function out = data_simulation(obj, freq)
% Generate the multi-static response (MSR) matrix of the Helmholtz equation (the perturbed field
% u-U) by solving the linear system (5.11) of [1]

for f=1:length(freq)
    % Compute the MSR matrix by evaluating the single layer potential
    MSR = zeros(obj.cfg.Ns_total, obj.cfg.Nr);
    
    % Resolution of the integral equation 5.11
    sol = PDE.Helmholtz_R2.solve_forward(obj.D{1}, freq(f), obj.cfg.all_src, obj.pmeb, obj.pmtt, obj.pmeb_bg, ...
        obj.pmtt_bg);
    vphi = sol(1:size(sol,1)/2, :);
    vpsi = sol(size(sol,1)/2+1:end, :);

    k0 = tools.Helmholtz.wavenb(freq(f), obj.pmeb_bg, obj.pmtt_bg);

    for s=1:obj.cfg.Ns_total
        rcv = obj.cfg.rcv(s); % receivers corresponding to the s-th source
        toto = ops.SingleLayer_H.eval(k0, obj.D{1}, vpsi(:,s), rcv);
        MSR(s,:) = transpose(toto);
    end
    
    out.MSR{f} = MSR;
    
    % Take the transpose so that each row corresponds to a source
    out.vphi{f} = transpose(vphi);
    out.vpsi{f} = transpose(vpsi);
end

out.freq = freq; % save the frequency list
end
