function [CGPTt, dt, Tmax, waveform] = theoretical_CGPT_time_recon(D, cnd, pmtt, ord, Tmax0, Ntime, scl, ovrspl)
% Compute the time-dependent CGPT matrix N = h*M (*:convolution) with h a
% waveform. The computation is done in the time domain by reconstruction.
% This is much faster compared to theoretical_CGPT_time
%
% Inputs:
% -D, cnd, pmtt: inclusions and their conductivity and permittivity constants.
% -ord: order of CGPT matrix
% -Tmax0: time duration of the pulse signal h (at the scale 1)
% -Ntime: number of resampling points equally distributed on [0, Tmax]. No resampling (only possible truncation) if
% Ntime=0 (default).
% -scl: scale of the pulse signal
% -overspl: oversampling factor 
%
% Outputs:
% -CGPTt: time-dependent CGPT matrices on [0, Tmax]
% -dt: time-step of CGPTt
% -Tmax: duration of the pulse signal at the given scale
% -waveform: the waveform of the pulse

Ntime0 = Ntime;

if nargin < 8
	ovrspl = 1;
end
ovrspl = max(1,floor(ovrspl));
Ntime = Ntime0 * ovrspl;

% The measurement circle must be sufficiently far from the object, in order
% to guarantee the decay of the truncation error
mradius = 6*(D.diameter/2+norm(D.center_of_mass));

K = max(20, 2*ord+1); % K (K>=5) is minimum number of sources for the reconstruction of the first order GPT
cfg = acq.Coincided([0,0]', mradius, K, [1, 2*pi, 2*pi], false, [1,-1]);

[waveform, dt, Tmax, ~] = tools.make_pulse(Tmax0, Ntime, scl);

P = PDE.PulseImaging_R2(D, cnd, pmtt, waveform, dt, cfg);

data = P.data_simulation();

% out = P.reconstruct_CGPT_analytic(data.MSR, ord);
out = P.reconstruct_CGPT_analytic(data.MSR, floor((K-1)/2));

CGPTt0 = zeros(2*ord, 2*ord, Ntime);

for t=1:Ntime
	CGPTt0(:,:,t) = out.CGPT{t}(1:2*ord, 1:2*ord);
end

CGPTt = CGPTt0(:,:,1:ovrspl:end);
end