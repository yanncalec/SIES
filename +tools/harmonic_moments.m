function mu = harmonic_moments(D, normal, Sigma, m, n)
%

% [COM, mass] = mass_center(D, normal, Sigma);
% xc = real(COM); yc = imag(COM);

mass = (sum(D(1,:).*normal(1,:).*Sigma) + sum(D(2,:).*normal(2,:).*Sigma))/2;
Cx = sum(1/2 * (D(1,:).^2) .* normal(1,:).*Sigma)/mass;
Cy = sum(1/2 * (D(2,:).^2) .* normal(2,:).*Sigma)/mass;

mu = zeros(m+1,n+1);
mu(1,1) = mass;
mu(2,1) = 0; mu(1,2) = 0;

for s=0:m
    for t=0:n
        if s>1 || t>1 || (s==1 && t==1)
            S1 = (D(1,:)-Cx).^(s+1) .* (D(2,:)-Cy).^t .* normal(1,:) / (s+1);
            S2 = (D(1,:)-Cx).^s .* (D(2,:)-Cy).^(t+1) .* normal(2,:) / (t+1);
        
            mu(s+1,t+1) = 1/2 * sum((S1 + S2).*Sigma);
        end
    end
end

