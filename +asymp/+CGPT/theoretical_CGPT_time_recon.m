function [CGPTt, dt, Tmax, waveform] = theoretical_CGPT_time_recon(D, cnd, pmtt, ord, Tmax0, Ntime, scl)
% Compute the time-dependent CGPT matrix N = h*M (*:convolution) with h a
% waveform. The computation is done in the time domain by reconstruction.
% This is much faster compared to theoretical_CGPT_time
%
% Inputs:
% D, cnd, pmtt: inclusions and their conductivity and permittivity
% constants.
% ord: order of CGPT matrix
% Tmax0: time duration of the pulse signal h (at the scale 1)
% Ntime: number of resampling points equally distributed on [0, Tmax]. No resampling (only possible truncation) if
% Ntime=0 (default).
% scl: scale of the pulse signal
%
% Outputs:
% CGPTt: time-dependent CGPT matrices on [0, Tmax]
% dt: time-step of CGPTt
% Tmax: duration of the pulse signal at the given scale

mradius = 2*(D.diameter/2+norm(D.center_of_mass));
K = 2*ord+1;
cfg = acq.Coincided([0,0]', mradius, K, [1, 2*pi, 2*pi], false, [1,-1]);

[waveform, dt, Tmax, ~] = tools.make_pulse(Tmax0, Ntime, scl);

P = PDE.PulseImaging_R2(D, cnd, pmtt, waveform, dt, cfg);

data = P.data_simulation();

out = P.reconstruct_CGPT_analytic(data.MSR, floor((K-1)/2));

CGPTt = zeros(2*ord, 2*ord, Ntime);

for t=1:Ntime
    CGPTt(:,:,t) = out.CGPT{t}(1:2*ord, 1:2*ord);
end

end