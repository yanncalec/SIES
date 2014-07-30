function [CGPTt, dt, CGPTf] = theoretical_CGPT_time(D, cnd, pmtt, ord, H, df, zp)
% Compute the time-dependent CGPT matrix N = h*M (*:convolution) with h a
% waveform. The computation is done in the frequency domain. 
%
% Inputs:
% D, cnd, pmtt: inclusions and their conductivity and permittivity
% constants.
% ord: order of CGPT matrix
% H: Fourier transform of the waveform h on [0, Fmax] with Fmax being the bandwidth of H. Note that the pulse h is a
% real function, hence H(-w) = conj(H(w)). We recall the convention for the Fourier transform:
%       H(w) = \int h(t) exp(-2*pi*1i*t*w) dt
% df: frequency sampling step for H
% zp: length of zero-padded H, 2^12 by default
%
% Outputs:
% CGPTt: time-dependent CGPT matrices on [0, Tmax] for some Tmax determined
% from H, df, and zp.
% dt: time-step of CGPTt
% CGPTf: CGPT matrix in the frequency domain of the band [-Fmax, Fmax] 

if nargin < 7
    zp = 2^12;
end

% Verify the maximum frequency so that |H(w) Mf(w)| <= eps
NH = length(H); % original length
M0 = asymp.CGPT.theoretical_CGPT(D, asymp.CGPT.lambda(cnd, pmtt, df*(NH-1)), ord);
toto = max(max(abs(H(end) * M0)));

if toto > 1e-5
    warning(['Enlarge the bandwidth of the waveform! The current residual is ', num2str(toto)]);
end

% H is the half (on [0, Fmax]) of FFT of h which has an even length, so the last term of H must be
% real. Add 0 if it is not the case. Zero-padding of H to a length 2^12, this increases the
% regularity of the result.
H = tools.zeropadding(H, zp, 'tail');
Nfreq = length(H); % Length in the time domain will be Ntime = 2*Nfreq
dt = 1/(df * Nfreq * 2); % time-step of CGPTt equals to 1/(2*Fmax)

% Compute Nf = H*Mf (*: product)
[nr, nc] = size(M0);
CGPTf = zeros(nr, nc, 2*Nfreq);

% If h is real and has even length N, then H(n)=H(N-n)^*, and precisely,
% H = [H(0), H(1)... H(N/2-1), H(N/2), H^*(N/2-1),... H^*(1)].
% For matlab, this becomes: H(n)=H(N-n+2), and
% H = [H(1), H(2)... H(N/2), H(N/2+1), H^*(N/2),... H^*(2)].

% So we compute only for the frequency in [0, Fmax]:
KsdS = asymp.CGPT.make_block_matrix(D);
if max(abs(pmtt)) == 0 % In this case the CGPT matrix is independent of the frequency
    M = asymp.CGPT.theoretical_CGPT_fast(D, KsdS, asymp.CGPT.lambda(cnd, pmtt, 0), ord);
    
    for f=1:NH % stop at NH since H=0 beyond NH (H zero-padded)
        CGPTf(:,:, f) = H(f) * M;
    end
else
    for f=1:NH % stop at NH since H=0 beyond NH (H zero-padded)
        CGPTf(:,:, f) = H(f) * asymp.CGPT.theoretical_CGPT_fast(D, KsdS, asymp.CGPT.lambda(cnd, pmtt, (f-1)*df), ord);
    end
end
%
for f=Nfreq+1 : 2*Nfreq % Copy to get the other part
    CGPTf(:,:,f) = conj(CGPTf(:,:,2*Nfreq-f+2));
end

% Take inverse Fourier transform to get h*M (*: convolution).
%
% No fftshift before ifft since CGPTf is alreaday shifted by
% construction. No fftshift after ifft since the signal is causal.

% ifft along the time axis
CGPTt = real(ifft(CGPTf, [], 3)) * (2*Nfreq) * df; % 2*Nfreq is due to the Matlab FFT convention, and df is due to the approximation of the true Fourier integral



