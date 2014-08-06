function [CGPTt, dt] = CGPT_time_truncation(CGPTt0, dt0, Tmax, Ntime)
% Truncate the 3D time dependent CGPT matrix on a time interval [0, Tmax]
% and resample it with Nt points
% Inputs:
% CGPTt0: time dependent CGPT 3D matrix (the third dimension is the time) 
% dt0: the time step in the CGPT
% Tmax: new time duration 
% Ntime: number of resampling points equally distributed on [0, Tmax]. No resampling (only possible truncation) if
% Ntime=0 (default).
% Output:
% CGPTt: new CGPT 3D matrix

if nargin < 4
    Ntime = 0;
end

Ntime0 = size(CGPTt0, 3);
Tmax0 = Ntime0 * dt0;

if Tmax > Tmax0
    error('Tmax must be smaller than the original time duration!');
end

T0 = (0:Ntime0-1) * dt0;

if Ntime <= 0 % no resampling, just truncation
    dt = dt0;
    Ntime = floor(Tmax/dt); 
    CGPTt = CGPTt0(:,:,1:Ntime);
else    
    dt = Tmax / Ntime;
    T1 = (0:Ntime-1) * dt;
    CGPTt = time_interp(CGPTt0, T0, T1);
end

function X1 = time_interp(X0, T0, T1)
% T0 and X0 must have the same length: length(T0)==size(X0,3)

[nr, nc, nt] = size(X0);

X1 = zeros(nr, nc, length(T1));

for rr=1:nr
    for cc = 1:nc
        X1(rr,cc,:) = interp1(T0, squeeze(X0(rr,cc,:)), T1, 'spline');
    end
end