function S = ShapeDescriptorFFP(W)
% Invariant shape descriptor computed from far field pattern
% Inputs:
% W: far field pattern (square) matrix. W can be 3d with frequency as the third dimension
%
% Outputs:
% S: rotation and translation invariant shape descriptor

nd = ndims(W);

if nd==2
    S = ShapeDescriptorFFP_onefreq(W);
    
elseif nd==3
    S = zeros(size(W));
    
    for n=1:size(W,3)
        S(:,:,n) = ShapeDescriptorFFP_onefreq(W(:,:,n));
    end
else
    error('The dimension of the input W must be 2 or 3!');
end

function S = ShapeDescriptorFFP_onefreq(W)
% Invariant shape descriptor computed from far field pattern

[Nv1, Nv2] = size(W);

if (Nv1~=Nv2)
    error('Input matrix must be square!');
end
Nv = Nv1;

G0 = abs(W);
G1 = fliplr(flipud(G0)); % Change of variable x -> -x;

S = real(ifft2(fft2(G0) .* fft2(G1))) * (2*pi/Nv)^2;

