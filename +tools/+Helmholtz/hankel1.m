function  y =hankel1(n, x)
% Hankel function of the first kind
% H^(1)_n(x)
% n is the order

y = besselj(n,x) + 1i*bessely(n,x);
