function [U] = cylind_wave(k, m, X)
% [U] = cylind_wave(k, m, X)
% Evaluate the cylindrical wave J_m(k*|X|)e^{im*theta_X} of given order m at given positions X.
% INPUTS:
% k: wave number, a scalar
% m: index of cylindrical wave 
% X: position of evaluation, 2-by-N
% OUTPUT
% U: a N-by-1 vector

    if size(X,1)~=2
        error('The inputs X must have two rows');
    end
    
    Z = X(1,:) + 1i*X(2,:);
    
    Jm = besselj(m, k*abs(Z));
    U = Jm .* exp(1i * m*angle(Z)); % cos(m*Angle) +1i*sin(m*Angle);
    U = U(:);
end
