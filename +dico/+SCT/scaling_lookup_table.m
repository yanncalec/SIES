function [err, scl, idx] = scaling_lookup_table(S0, frange, S, sfrange)
% Given the shape-descriptor S, find the scaling factor by looking up the reference table of
% shape-descriptor S0.
%
% INPUTS:
% -S0: frequency dependent shape-descriptor (inv. to translation and rotation)
% -frange: frequency range of S0 [f1, f2]
% -S: shape-descriptor extracted from multi-frequency data
% -sfrange: scanning frequency range of S (used for data acquisition) [sf1, sf2]
% OUTPUTS:
% -err: euclidean distance between S and S0 for different scaling factor in the table
% -scl: values of scaling factors in the table
% -idx: index of the best scaling factor minimizing err

% range of acceptable scaling factor [s1,s2]
sclrange = frange ./ sfrange;

Nf = length(S0);

% resolution of frequency in S0
dfreq = (frange(2)-frange(1))/Nf;

% resolution of scaling factor
dscl = dfreq / sfrange(2) / 4;
% dscl = dfreq / sfrange(1);

Ns = ceil((sclrange(2)-sclrange(1))/dscl);
err = zeros(1,Ns);

% table of values of scaling factor
scl = linspace(sclrange(1),sclrange(2),Ns);

for n=1:Ns
	% index in the table
	idx = floor((scl(n)*sfrange - frange(1))/dfreq);
	idx = [max(idx(1),1), min(idx(2), Nf)];
	
	err(n) = tools.compare_sounds(S0(idx(1):idx(2)), S);
end

[~, idx] = min(err);
idx = idx(1);

% %%
% Ns = floor((sclrange(2)-sclrange(1))/dscl);
% err = zeros(1,Ns+1);

% % table of values of scaling factor
% scl = sclrange(1) + (0:Ns)*dscl;

% for n=1:(Ns+1)
%     % index in the table
%     idx = floor(scl(n)*sfrange/dfreq);
%     idx = [max(idx(1),1), min(idx(2), Nf)];

%     err(n) = tools.compare_sounds(S0(idx(1):idx(2)), S);
% end

% [~, idx] = min(err);
% idx=idx(1);
