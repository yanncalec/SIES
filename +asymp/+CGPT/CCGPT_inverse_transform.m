function [Z1, Z2] = CCGPT_inverse_transform(Y1, Y2, T0, S0, Phi0)
% [Z1, Z2] = CCGPT_inverse_transform(Y1, Y2, T0, S0, Phi0)
%
% On the CCGPT data Y1, Y2, apply first the translation -T0, then the
% scaling 1/S0, and the rotation -Phi0. The outputs Z1, Z2 are the
% transformed data.
%

if nargin<5
    Phi0 = 0;
end

ord = size(Y1, 1);
R = repmat((1:ord)', 1, ord);
Comb = zeros(ord);

for n=1:ord
    for m=1:n
        Comb(n,m) = nchoosek(n,m);
    end
end

if length(T0)>1 % if a vector is passed
    T0 = (T0(1) + T0(2) * 1j);
end

% Ct * Comb .* ((T0).^(R - R')) is identity
Ct = Comb .* ((-T0).^(R - R'));  
Ct(isnan(Ct)) = 0; Ct(isinf(Ct)) = 0;

Gy = diag((1/S0*exp(-1j*Phi0)).^(1:ord));

Z1 = Gy * Ct * Y1 * Ct.' * Gy;
Z2 = conj(Gy) * conj(Ct) * Y2 * Ct.' * Gy;

% $$$ Am = diag(exp(-j*Phi0).^(1:ord));
% $$$ Sm = diag((1/S0).^(1:ord));
% $$$ 
% $$$ Z1 = Sm * Am * Ct * Y1 * Ct.' * Am * Sm;
% $$$ Z2 = Sm * conj(Am) * conj(Ct) * Y2 * Ct.' * Am * Sm;

