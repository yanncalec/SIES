function [ Y1 ] = P1_interp(Y, m, idx)
% Interpolation with P1 elements
% Y: coefficients of P1 elements
% m: sampling step of P1 elements
% idx: index of points where the interpolated value is desired

N = length(Y) * m;
Psi = tools.BEM.P1_basis(N, m);
V = sum(diag(Y) * Psi, 1);
if ~isempty(idx)
    Y1 = V(idx);
else
    Y1 = V;
end

end