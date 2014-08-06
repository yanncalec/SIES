function SD = ShapeDescriptor_PT_time(CGPT, Scl, dt, Tmax1)
% Inputs:
% CGPT: a cell of 3D time-dependent CGPT matrix, CGPT{i} is the CGPT at
% the i-th scale.
% Scl: scaling parameter for each scale
% dt: time step of CGPT, optional
% Tmax1: time duration of CGPT at the 1st scale, optional
% Output:
% SDt: shape descriptor

if nargin < 4
    Tmax1 = 1;
end

if nargin < 3
    dt = 1;
end

if ~iscell(CGPT) % transform to a cell
    CGPT = {CGPT};
end

scl = length(CGPT); % total number of scales
Ntime = size(CGPT{1}, 3);

Ft = zeros(Ntime, scl);

for s = 1:scl
    for t=1:Ntime
        Ft(t, s) = norm(squeeze(CGPT{s}(1:2,1:2,t)), 'fro');
    end
    
    % renormalization, since the pulse waveform at the scale s is h_s(t) = Scl(s)*h(Scl(s)*t)
    Ft(:, s) = Ft(:, s) / Scl(s);    
end

% Invariant to dilation:
% Renormalization using the first scale information 
cst = sum(Ft(:, 1)) * dt / Tmax1; 
SD = reshape(Ft/cst, 1, []);

end

