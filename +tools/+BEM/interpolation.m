function [ Y1 ] = interpolation(Psi, Y, idx)
% [ Y1 ] = interpolation(Psi, Y, idx)
% Interpolation using a basis Psi, with coefficient vector Y, return the value
% evaluated at positions of index idx. Psi is a N-by-M matrix with M the dimension
% of Y, and each column of Psi contains the value of one basis function psi
% evaluated at a fine lattice of N points.
%
    if nargin < 3
        idx = [];
    end
    
    V = Psi * reshape(Y, size(Psi, 2), []);
    if ~isempty(idx)
        Y1 = V(idx,:);
    else
        Y1 = V;
    end
end

