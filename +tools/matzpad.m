function V = matzpad(X, dim)
% V = matzpad(X, dim)
% Extend a 2D matrix to dimension dim by zero-padding

[m, n] = size(X);
if dim(1)< m || dim(2) < n
    error('Dimension error');
end

V = zeros(dim);
r0 = floor(dim(1)/2)+1;
c0 = floor(dim(2)/2)+1;

r1 = r0 - floor(m/2);
c1 = c0 - floor(n/2);

V(r1:r1+m-1, c1:c1+n-1) = X;

