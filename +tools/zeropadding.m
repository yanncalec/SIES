function [Xz] = zeropadding(X, n)
% Zero pad a vector X to a length of power of 2^n.

Nx = length(X);
Nz = 2^n;

if Nz < Nx
    Xz = X;
    %error('Invalid factor n!');
else
    Xz = zeros(Nz, 1); Xz(1:Nx) = X;
end




