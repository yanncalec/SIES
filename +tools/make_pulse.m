function [waveform, dt, Tmax, freqform, df, Fmax, extrema] = make_pulse(Tmax0, Ntime, scl, L2nrm)
% Make gaussian pulse h (derivative of gaussian) at a given scale and its Fourier transform H.
%
% Inputs:
% Tmax0: time duration of the initial pulse signal h (at the scale 1)
% Ntime: (minimal) number of time steps on [0, Tmax]
% scl: scaling parameter. The output pulse is h_s = scl*h(scl*t)
% L2nrm: if true then the pulse is normalized to keep the constant L^2
% energy, optional.
%
% Outputs:
% waveform: discrete values of the pulse in [0, Tmax]
% dt: dt=Tmax/Ntime, constant independent of ord, scl
% Tmax: real duration of the output waveform
% freqform: H(w) for discrete values of w in [0, Fmax], Fmax is the frequency bandwidth so that abs(H(Fmax))<1e-8
% df: frequency step of freqform, constant independent of ord, scl
% Fmax: half band width of the output waveform (larger than 1e-8)
% extrema: index n where n*dt is a extrema of h

% Convention for Fourier transform:
%       f^(w) := \int f(t) e^(-2pi*t*w) dt.
% Fourier transform of exp(-a*x^2) is
%   \sqrt(pi/a) * exp(-pi^2*w^2/a)
% Fourier transform of x^d*exp(-pi*x^2) is
%   (exp(-pi*w^2))^(d) / (-2*pi*i)^d

if nargin < 4
    L2nrm = false;
end

if nargin < 3
    scl = 1;
end

if scl<0
    error('Scaling must be positive!');
end

if nargin < 2
    Ntime = 2^10;
end

if nargin < 1
    Tmax0 = 5; % time duration of the pulse h(x) = (exp(-pi*(x-T0)^2))^{(3)} (to be symmetric about 0)
end    

if Tmax0 < 5
    warning('Truncation (Tmax) for the initial waveform is too small!');
end

syms x real positive
% syms x real 
% assume(x>=0);

T0 = 2.5; % -T0 is the starting time (ie, h(0) is considered as 0)
ord = 3; % order of the derivative

Tmax = Tmax0 / scl; % duration at the scale 1 is about 5

f = exp(-pi * (x-T0)^2); 
h = simplify(diff(f, x, ord)); % waveform h is the d-th derivative of f
% f = x^d * exp(-pi * x^2);

wavefunc = scl * subs(h, x, scl*x); % hj(x) = scl * h(scl * x), L^1 normalisation
H = (2*pi*1i*x)^ord * exp(-2*pi*1i*T0*x) * exp(-pi*x^2); % Fourier transform of h
% H = simplify((2*pi*1i*x)^ord / (-2*pi*1i)^d * diff(exp(-pi*x^2), d)); % Fourier transform of h
freqfunc = subs(H, x, x/scl);

% Estimate the essential bandwidth of H
% freqfunc0 = (2*pi*x)^ord * exp(-pi*x^2);
S = solve(abs(freqfunc) - 1e-8, 'Real', true); % keep only real solution
Fmax = max(double(S)); % Half bandwidth of H = h^.

% We truncate h on [0, Tmax] (and H on [0, Fmax]), and sample h (and H)
% using Ntime (and Nfreq = Ntime) points, such that Ntime >> 2*Fmax*Tmax.
Ntime = max(Ntime, 2*Tmax*Fmax); 

% Time domain evaluation
dt = Tmax/Ntime; % time step
hj = inline(wavefunc);
waveform = hj((0:Ntime-1)*dt);
% Or equivalently (but it is much slower than the inline function)
% waveform = double(subs(wavefunc, x, (0:Ntime-1)*dt)); 

% Frequency domain evaluation
Nfreq = Ntime;
df = Fmax/Nfreq;
% df = min(1/2/Tmax, 10^-2) / 2/Tmax; 
% Nfreq = ceil(Fmax/df);
Hj = inline(freqfunc);
freqform = Hj((0:Nfreq-1)*df);

dh = simplify(diff(h,x,1));
E = solve(dh, 'Real', true)/scl;
extrema0 = double(sort(ceil(E/dt)));
extrema = extrema0(extrema0<Ntime);

% For L^2 normalization
if L2nrm
    waveform = waveform / sqrt(scl);
    freqform = freqform / sqrt(scl);
end

% In case that the analytical form of H is not known, we can also compute by fft as follows:
%
% In order to get a high precision approximation, we increase
% Tmax (the frequency step is 1/Tmax) and Ntime (which reduces
% the time-step Tmax/(Ntime-1)).
%
% Ntime0 = max(Ntime, 2^13); Tmax0 = dt*Ntime0; % new Tmax and Ntime
% df = 1/Tmax0; % frequency step
%
% waveform0 = hj((0:Ntime0-1)*dt);
% freqform0 = dt * fft(waveform0); % Fourier transform (approximation by FFT)
%
% Keep the low frequency of [0, Fmax]. The frequency beyond
% Fmax is close to zero. We remove them to accelerate the
% computation in theoretical_CGPT_time (computation of CGPT for
% each frequency).
%
% Number of entries to keep, must be even
% Nc = min(Ntime0/2, ceil(2*Fmax/df/2))*2; % 2*Fmax/df
% freqform = freqform0(1:Nc/2); % the term Nc/2+1 is real
end
