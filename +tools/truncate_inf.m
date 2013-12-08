function Y = truncate_inf(X)
idx1 = isinf(X);
idx2 = isfinite(X);
sgn = sign(X(idx1));
fmax = max(abs(X(idx2)));
Y = X;
Y(idx1) = fmax*sgn;
