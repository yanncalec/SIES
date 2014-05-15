function [CGPTt] = CGPT_time_truncation(CGPTt0, dt, Tmax)

[nr, nc, Ntime0] = size(CGPTt0);

% We keep finally [0, Tmax]
Ntime = ceil(Tmax/dt); % 2*Tmax/Tmax0*(Nfreq0-1) = Tmax/Tmax0*Ntime0
CGPTt = zeros(nr, nc, Ntime);

if Ntime > Ntime0 % add zeros in this case
    CGPTt(:,:,1:Ntime0) = CGPTt0;
else
    CGPTt = CGPTt0(:,:,1:Ntime);
end
