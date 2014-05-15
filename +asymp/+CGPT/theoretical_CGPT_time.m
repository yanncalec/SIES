function [CGPTt, dt, CGPTf] = theoretical_CGPT_time(D, cnd, pmtt, ord, H, df)
% Compute the time-dependent CGPT matrix N = h*M (*:convolution) with h a
% waveform. The computation is done in the frequency domain. 
%
% Inputs:
% D, cnd, pmtt: inclusions and their conductivity and permittivity
% constants.
% ord: order of CGPT matrix
% H: Fourier transform of the waveform h on [0, Fmax]. Note that h is a
% real function, hence H(-w) = conj(H(w)).
% We recall the convention for the Fourier transform:
%       H(w) = \int h(t) exp(-2*pi*1i*t*w) dt
% df: frequency step of H
% Tmax: stopping time
%
% Outputs:
% CGPTt: time-dependent CGPT matrices on [0, Tmax] for some Tmax determined
% from H and df, Tmax = dt*length(H)
% dt: time-step for Nt

% Verify the maximum frequency so that |H(w) Mf(w)| <= eps
NH = length(H);
M0 = asymp.CGPT.theoretical_CGPT(D, asymp.CGPT.lambda(cnd, pmtt, df*(NH-1)), ord);
toto = max(max(abs(H(end) * M0)));

if toto > 1e-7
    Warning('Enlarge the bandwidth of the waveform!');
end

% H is the half (on [0, Fmax]) of FFT of h which has an even length, so the last term of H must be
% real. Add 0 if it is not the case. Zero-padding of H to a length 2^12, this increases the
% regularity of the final result.
H = tools.zeropadding(H, 12); 
Nfreq0 = length(H);

Ntime0 = 2*Nfreq0; % Length in the time domain

% Fmax: bandwidth of H. It must be sufficiently large to reduce the
% numerical error (by keeping H(w)Mf(w) small)
Fmax = df * Nfreq0; % Nfreq0/Tmax0, with Tmax0 defined in the function make_pulse 
dt = 1/2/Fmax; % time-step

% Compute Nf = H*Mf (*: product)
[nr, nc] = size(M0);

CGPTf = zeros(nr, nc, Ntime0);

% If h is real and has even length N, then H(n)=H(N-n)^*, and precisely,
% H = [H(0), H(1)... H(N/2-1), H(N/2), H^*(N/2-1),... H^*(1)].
% For matlab, this becomes: H(n)=H(N-n+2), and
% H = [H(1), H(2)... H(N/2), H(N/2+1), H^*(N/2),... H^*(2)].

% So we compute only for the frequency in [0, Fmax]:
for f=1:NH % stop at NH since H=0 beyond NH (H zero-padded)
    freq = (f-1) * df;
    toto = asymp.CGPT.theoretical_CGPT(D, asymp.CGPT.lambda(cnd, pmtt, freq), ord);
    CGPTf(:,:, f) = toto * H(f);
end
%
for f=Nfreq0+1 : 2*Nfreq0 % Copy to get the other part
    CGPTf(:,:,f) = conj(CGPTf(:,:,2*Nfreq0-f+2));
end

% Take inverse Fourier transform to get h*M (*: convolution).
%
% We compute f(m*dt) for m=0..Nfreq-1 with dt=1/2*Fmax, via FFT:
% f(m*dt) \simeq df * \sum_{n=-Nfreq/2}^{Nfreq/2-1} F(n*df) *
% exp(2*pi*1i*df*n*m*df*dt), note that df*dt = 1/(Nfreq-1)
%            =    2*Fmax * 1/2(Nfreq-1) \sum_{n=0}^{Nfreq-1} F(n*df)
%            * exp(2*pi*1i*n*m/(Nfreq-1))
%            =    2*Fmax * ifft(fftshift(F))
% There is no fftshift after ifft since the index m varies in 0..Nfreq-1
% (Nt is supported on R_+).

CGPTt = zeros(nr, nc, Ntime0);
%
for rr=1:nr
    for cc=1:nc
        y = squeeze(CGPTf(rr,cc,:)); % No fftshift after y since CGPTf is alreaday shifted by construction
        CGPTt(rr,cc,:) = real(ifft(y)) * 2 * Fmax; % Remind the definition of ifft in matlab!
    end
end
