function [S, G] = ShapeDescriptorSCT(W, Nv, bdwidth)
% [S, G] = ShapeDescriptorSCT(W, Nv, bdwidth)
% Invariant shape descriptor based on Scattering coefficients
% Inputs:
% W: SCT (square) matrix. W can be 3d with frequency as the third dimension
% Nv: sampling points for the square [0,2pi]^2, Nv=256 by default
% bdwidth: if >0 then compute with band-diagonal (of width bdwidth) far feild pattern
%
% Outputs:
% S: rotation and translation invariant shape descriptor
% G: far field pattern

if nargin<3 || bdwidth == 0
    bdmask=[];
else
    bdmask = tools.bandiag_mask(Nv, bdwidth);
end

if nargin<2
    Nv =2^8;
end

nd = ndims(W);

if nd==2
    [S, G] = ShapeDescriptorSCT_onefreq(W, Nv, bdmask);

elseif nd==3
    [~,~,Nf] = size(W);

    S = zeros(Nv,Nv,Nf);
    G = zeros(Nv,Nv,Nf);
    
    for n=1:Nf
        [S(:,:,n), G(:,:,n)] = ShapeDescriptorSCT_onefreq(W(:,:,n), Nv, bdmask);
    end
else
    error('The dimension of the input W must be 2 or 3!');
end

function [S, G] = ShapeDescriptorSCT_onefreq(W, Nv, bdmask)
% Invariant shape descriptor based on Scattering coefficients
% Inputs:
% W: SCT (square) matrix
% Nv: sampling points for the square [0,2pi]^2, Nv=256 by default
% bdmask: band-diagonal mask, not used by default
%
% Outputs:
% S: rotation and translation invariant shape descriptor
% G: far field pattern

G = dico.Helmholtz.farfieldpattern(W,Nv);

if nargin == 3 || ~isempty(bdmask)
    G = G .* bdmask;
end

G0 = abs(G);

G1 = fliplr(flipud(G0)); % Change of variable x -> -x                        
S = real(ifft2(fft2(G0) .* fft2(G1))) * (2*pi/Nv)^2;

% % Strangely, conv2 is extremely slow and inaccurate
% tic
% S = conv2(G0, G1, 'same') * (2*pi/Nv)^2;
% toc
