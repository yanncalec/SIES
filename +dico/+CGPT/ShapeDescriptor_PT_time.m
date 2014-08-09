function SDt = ShapeDescriptor_PT_time(CGPT, Scl)
% Inputs:
% CGPT: a cell of 3D time-dependent CGPT matrix, CGPT{i} is the CGPT at
% the i-th scale.
% Scl: scaling parameter for each scale
% Output:
% SDt: shape descriptor of size (Ntime X length(Scl))

if ~iscell(CGPT) % transform to a cell
    CGPT = {CGPT};
end

scl = length(CGPT); % total number of scales
Ntime = size(CGPT{1}, 3);

SD = zeros(Ntime, scl, 2);

for s = 1:scl
    for t=1:Ntime
        SD(t, s, :) = svd(squeeze(CGPT{s}(1:2,1:2,t)));
    end
    
    % renormalization, since the pulse waveform at the scale s is h_s(t) = Scl(s)*h(Scl(s)*t)
    SD(:, s, :) = SD(:, s, :) / Scl(s);
end

% Invariant to dilation:
% Renormalization using the first scale information 
% Old version:
% toto = squeeze(SD(:, 1, :));
% cst = mean(toto, 1);
% 
% SDt = SD / sum(cst);
% % SDt = SD / cst(2); % this can also work

SDt = zeros(Ntime, scl);
for s = 1:scl
    SDt(:, s) = SD(:,s,1) / mean(SD(:,1,2));
end

% Other possible ways to construct invariants: 
% SD1 = S(:,:,1)/cst(1); 
% SD2 = S(:,:,2)/cst(2);
%
% Or:
%
% SD1 = (S(:,:,1)+S(:,:,2))/sum(cst);
% SD2 = SD1;
%
% But these don't seem to work: confusion of ellipse with flower
end

