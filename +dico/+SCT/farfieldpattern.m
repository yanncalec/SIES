function [G] = farfieldpattern(W, Nv)
% Inputs:
% -W: SCT (square) matrix
% -Nv: sampling points for the square [0,2pi]^2, Nv=256 by default
% Outputs:
% -G: far field pattern

if nargin < 2
	Nv = 256;
end

[M,N] = size(W);

if M~=N
	error('Input SCT matrix must be square!');
end

% Compute the Fourier series:
W1 = tools.matzpad(W, [Nv, Nv]); % zero-padding

% F = fftshift(fft2(W1)); % Fourier series
F = fft2(W1); % Fourier series

G0 = fliplr(F);

tt=floor(Nv/4);
G=circshift(G0, [tt, tt]);;
