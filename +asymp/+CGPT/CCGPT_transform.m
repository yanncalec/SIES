function [Z1, Z2] = CCGPT_transform(Y1, Y2, T0, S0, Phi0)
% [Z1, Z2] = CCGPT_transform(Y1, Y2, T0, S0, Phi0)
%
% On the CCGPT Y1, Y2, apply first the rotation Phi0, the scaling S0,
% then translation T0. The outputs Z1, Z2 are the transformed CCGPT.
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

t0 = (T0(1) + T0(2) * 1j);
Ct = Comb .* (t0.^(R - R'));
Ct(isnan(Ct)) = 0; Ct(isinf(Ct)) = 0;

Gy = diag((S0*exp(1j*Phi0)).^(1:ord));

Z1 = Ct * Gy * Y1 * Gy * (Ct.');
Z2 = conj(Ct) * conj(Gy) * Y2 * Gy * (Ct.');
