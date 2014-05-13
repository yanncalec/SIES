function [Nt, timeStep] = theoretical_CGPT_time(D, cnd, pmtt, ord, H, fmax)
% Compute the time-dependent CGPT matrix N = h*M (*:convolution) with h a
% waveform. The computation is done in the frequency domain. 
%
% Inputs:
% D, cnd, pmtt: inclusions and their conductivity and permittivity
% constants.
% ord: order of CGPT matrix
% H: Fourier transform of the waveform h on [0, fmax]. Note that h is a
% real function, hence H(-w) = conj(H(w)). Recall the convention:
%   H(w) = \int h(t) exp(-2*pi*1i*t*w) dt
% fmax: bandwidth of H. fmax must be sufficiently large to reduce the
% numerical error (by keeping H(w)Mf(w) small)


% Verify the maximum frequency so that |H(w) Mf(w)| < 1e-8
epsilon = 1e-8;

lambda0 = asymp.CGPT.lamda(cnd, pmtt, fmax);
M0 = asymp.CGPT.theoretical_CGPT(D, lambda0, ord);
toto = max(abs(H(end) * M0));

if toto > epsilon || 1/fmax > 1e-2
    error('Enlarge the bandwidth of the waveform!');
end

% Get H on [-fmax, fmax]
Nfreq0 = length(H);
df = fmax / (Nfreq0-1); % frequency-step

freq = (1-Nfreq0 : Nfreq0-2) * df;
Nfreq = length(freq); % 2*(Nfreq0 - 1)

lambda = asymp.CGPT.lamda(cnd, pmtt, freq);
Hv = [reshape(conj(H(end:-1:2)), [], 1); H(1:end-1)];

% Compute Nf = H*Mf (*: product)
Mf = {};
[nr, nc] = size(M0);

Nf = zeros(nr, nc, Nfreq);

for f=1:Nfreq
    toto = asymp.CGPT.theoretical_CGPT(D, lambda(freq(f)), ord); 
    Nf(:,:, f) = toto * Hv(f);
end

% Take inverse Fourier transform to get N = h*M (*: convolution)
Nt = zeros(nr, nc, Nfreq);
% dt = 1/2/fmax; % time-step

% By our choice of dt and df, it holds dt*df = 1/(Nfreq-1)

% f(m*dt) \simeq df * \sum_{n=1-Nfreq0}^{Nfreq0-2} F(n*df) exp(2*pi*1i*df*n*m*dt)
%            =    2*fmax * 1/2(Nfreq0-1) \sum_{n=1}^{2(Nfreq0-1)} F(n*df)
%            exp(2*pi*1i*n*m/2(Nfreq0-1))
%            =    2*fmax * fftshift(ifft(fftshift(F)))
for rr=1:nr
    for cc=1:nc
        y = squeeze(Nf(rr,cc,:));
        Nt(rr,cc,:) = ifft(fftshift(y)) * 2 * fmax; 
        % There is no the second fftshift since Nt is supported on R_+
    end
end

timeStep = (0 : 2*(Nfreq0-1)-1)/2/fmax;




