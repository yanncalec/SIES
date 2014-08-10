function SDt = ShapeDescriptor_PT_time(CGPT, Scl, method)
% Inputs:
% CGPT: a cell of 3D time-dependent CGPT matrix, CGPT{i} is the CGPT at
% the i-th scale.
% Scl: scaling parameter for each scale
% Output:
% SDt: shape descriptor of size (Ntime X length(Scl))

if nargin < 3
    method = 1;
end

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

% Invariant to dilation: Renormalization using the first scale information. 
% The idea is to use the ratio of the two singular values (larger/smaller). 
% Mathematically however, one can prove that (larger/sqrt(smaller^2+larger^2)) is always
% well defined. Numerically these two methods are similar.
% The smaller singular value doesn't bring enough information: the
% rotaional symmetric objects can be very close.

if method == 1
    % Method 1: L1 norm
    cst = mean(squeeze(SD(:, 1, :)), 1);
    SDt = squeeze(SD(:,:,1) / sum(cst));
elseif method == 2
    % Method 1: L2 norm
    cst = mean(sqrt(squeeze(SD(:, 1, 1).^2 + SD(:, 1, 2).^2)));
    SDt = squeeze(SD(:,:,1) / cst);
elseif method == 3
    % Method 3: simple ratio ( this is more sensitive to noise )
    SDt = squeeze(SD(:,:,1) / mean(SD(:,1,2)));
elseif method == 4
    % Method 4: use two singular values
    cst = mean(squeeze(SD(:, 1, :)), 1);
    SDt1 = squeeze(SD(:,:,1) / sum(cst));
    % SDt1 = squeeze(SD(:,:,1) / mean(SD(:,1,2)));
    
    cst = mean(sqrt(squeeze(SD(:, 1, 1).^2 + SD(:, 1, 2).^2)));
    SDt2 = squeeze(SD(:,:,1) / cst);
    
    SDt = zeros(Ntime, scl, 2); SDt(:,:,1) = SDt1; SDt(:,:,2) = SDt2;
    % SDt = (SDt1 + SDt2)/2;
elseif method == 5
    % Method 5: use two singular values
    cst = mean(squeeze(SD(:, 1, :)), 1);
    SDt = SD / sum(cst);
elseif method == 6
    % Method 6: use two singular values with L2 renormalization
    cst = mean(sqrt(squeeze(SD(:, 1, 1).^2 + SD(:, 1, 2).^2)));
    SDt = SD / cst;
else
    error('Unknown method!');
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
