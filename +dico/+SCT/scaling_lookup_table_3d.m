function [err, scl] = scaling_lookup_table_3d(S0, frange, S, sfrange)
% Given the shape-descriptor S, find the scaling factor by looking up the reference table of
% shape-descriptor S0.
%
% INPUTS:
% -S0: pre-computed frequency dependent shape-descriptor (inv. to translation and rotation), a 3D
% array. The third dimension of S0 is the frenquency.
% -frange: frequency range of S0 [f1, f2]
% -S: shape-descriptor extracted from multi-frequency data, a 3D array.
% -sfrange: scanning frequency range of S (used for data acquisition) [sf1, sf2]
% OUTPUTS:
% err: euclidean distance between S and S0 for different scaling factor in the table
% scl: values of scaling factors in the table


% range of acceptable scaling factor [s1,s2]
sclrange = frange ./ sfrange

Nf = size(S0, 3);
% resolution of frequency in the table
dfreq = (frange(2)-frange(1))/Nf

% resolution of scaling factor
% dscl = dfreq / sfrange(2) 
dscl = dfreq / sfrange(1)

Ns = floor((sclrange(2)-sclrange(1))/dscl)
err = zeros(1,Ns+1);

for n=0:Ns
    % scaling factor to look up in the table
    scl = sclrange(1) + n*dscl; 

    % index in the table
    idx = floor(scl*sfrange / dfreq); 
    idx = [max(idx(1),1), min(idx(2), Nf)]

    V = S0(:,:,idx(1):idx(2));
    err(n+1) = tools.compare_movies(V, S);
end

[v, idx] = min(err);
scl = sclrange(1) + (idx-1)*dscl; 

% function compare_by_interpl(X,Y)

% X=X(:); Y=Y(:);
% N = length(X) - length(Y);
% Y = [Y; zeros(N,1)];
% fX = fft(X); fY = fft(Y);
