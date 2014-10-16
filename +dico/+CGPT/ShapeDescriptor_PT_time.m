function [SDt, SD] = ShapeDescriptor_PT_time(CGPT, method, Scl)
% Compute the time-dependent shape descriptor.
% Inputs:
% -CGPT: a cell of 3D time-dependent CGPT matrix, CGPT{i} is the CGPT at
% the i-th scale.
% -method: method for computing the shape descriptor, see the code
% -Scl: normalisation of the shape descriptor (to have the same numerical
% values across scales), optional. One must be careful in the convention 
% of normalization for the shape descriptor (eg, normalization by sqrt(Scl) 
% corresponds to the case L2nrm=true in the function tools.make_pulse).
%
% Output:
% -SDt: shape descriptor of size (Ntime X length(Scl))
% -SD: normalized singular values (based on which SDt is computed)

nbScl = length(CGPT); % total number of scales

if nargin < 3
    Scl = ones(1,nbScl);
end

if nargin < 2
    method = 2;
end

if ~iscell(CGPT) % transform to a cell
    CGPT = {CGPT};
end

Ntime = size(CGPT{1}, 3);

SD = zeros(Ntime, nbScl, 2);
FD = zeros(Ntime, nbScl);

for s = 1:nbScl
    for t=1:Ntime
        SD(t, s, :) = svd(squeeze(CGPT{s}(1:2,1:2,t))) / Scl(s);
        FD(t, s) = norm(CGPT{s}(1:2,1:2,t), 'fro')^2;
    end    
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
    cst = sqrt(mean(squeeze(SD(:, 1, 1).^2 + SD(:, 1, 2).^2)));
    SDt = squeeze(SD(:,:,1) / cst);
elseif method == 3
    % Method 3: simple ratio ( this is more sensitive to noise )
    SDt = squeeze(SD(:,:,1) / mean(SD(:,1,2)));
elseif method == 4
    % Method 4: use two singular values
    cst = mean(squeeze(SD(:, 1, :)), 1);
    SDt1 = squeeze(SD(:,:,1) / sum(cst));
    % SDt1 = squeeze(SD(:,:,1) / mean(SD(:,1,2)));
    
    cst = sqrt(mean(squeeze(SD(:, 1, 1).^2 + SD(:, 1, 2).^2)));
    SDt2 = squeeze(SD(:,:,1) / cst);
    
    SDt = zeros(Ntime, nbScl, 2); SDt(:,:,1) = SDt1; SDt(:,:,2) = SDt2;
    % SDt = (SDt1 + SDt2)/2;
elseif method == 5
    % Method 5: use two singular values
    cst = mean(squeeze(SD(:, 1, :)), 1);
    SDt = SD / sum(cst);
elseif method == 6
    % Method 6: use two singular values with L2 renormalization
    cst = sqrt(mean(squeeze(SD(:, 1, 1).^2 + SD(:, 1, 2).^2)));
    SDt = SD / cst;
elseif method == 7
    cst = mean(FD(:,1));
    SDt = sqrt(FD / cst);
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
