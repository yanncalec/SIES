function [U] = plane_wave(k, X, Xi)
% [U] = plane_wave(k, X, Xi)
% Evaluate the plane wave exp(k*i*<X,Y>) of given directions Xi at given positions X.
% INPUTS:
% k: wave number
% X: position of evaluation, 2-by-N
% Xi: direction of plane waves, 2-by-M with normalized columns
% OUTPUT
% U: a N-by-M matrix, with U(n,m)=exp(1i*k*<X(:,n), Xi(:,m)>)
    
    if size(X,1)~=2 || size(Xi,1)~=2
        error('The inputs X and Xi must have two rows!');
    end
    
    toto = tools.tensorprod(X(1,:), Xi(1,:)) + tools.tensorprod(X(2,:), Xi(2,:));
    U = exp(1i*k*toto);
end
