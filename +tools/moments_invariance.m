function phi = moments_invariance(mu)

[M,N] = size(mu);
if M<4 || N<4
    error('Moments of order 3 are needed!');
end

mu00 = mu(1,1);
mu = diag(mu00.^(-(1:M)/2)) * mu * diag(mu00.^(-(1:N)/2));

mu20 = mu(3,1);
mu02 = mu(1,3);
mu11 = mu(2,2);
mu30 = mu(4,1);
mu03 = mu(1,4);
mu12 = mu(2,3);
mu21 = mu(3,2);


phi(1) =  mu20 + mu02;
phi(2) = (mu20 - mu02)^2 + 4*mu11^2;
phi(3) = (mu30 - 3*mu12)^2 + (3*mu21 - mu03)^2;
phi(4) = (mu30 - mu12)^2 + (mu21 - mu03)^2;
