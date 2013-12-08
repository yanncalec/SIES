function [C, f, df, ddf] = poly_petal_interp(epsilon,Npetal,tau, theta)
% [C, f, df, ddf] = poly_petal_interp(epsilon,Npetal,tau, theta)
% Calculate a 6 order polynomial to simulate a dammaged petal
% Inputs:
% epsilon, Npetal: the parametres of make_flower
% tau: the percentage of dommage, eg, tau=0.1 means 90% of petal is
% good
% theta: the boundary parameterization of the first petal
%
% Outputs:
% f,df,ddf: the polynomial value, its first and second derivative
% on theta

w = 2*pi/Npetal;

%f(0)=f(w)=1-epsilon;
a = 1-epsilon;

%f(w/2) = tau*(1+epsilon);
%b = tau*(1+epsilon);
b = 1+epsilon*(2*tau-1);

%f'(0)=f'(w)=0
%f''(0)=f''(w)=Npetal^2*epsilon
c = Npetal^2*epsilon;

a0 = b;

% a2, a4, a6 are determined by solving a 3X3 linear system:

A = [2^(-6)*w^4, 2^(-4)*w^2, 1/4; 6*2^(-5)*w^4, 4*2^(-3)*w^2, 1; ...
     30*2^(-4)*w^4, 12*2^(-2)*w^2, 2];

Y = [(a-b)/w^2, 0, c]';

Coeff = inv(A)*Y;
C = [Coeff; a0];

a6 = Coeff(1); a4 = Coeff(2); a2 = Coeff(3); 

theta = theta-w/2;

f = a6*theta.^6 + a4*theta.^4 + a2*theta.^2 + a0 ;

df = 6*a6*theta.^5 + 4*a4*theta.^3 + 2*a2*theta ;

ddf = 30*a6*theta.^4 + 12*a4*theta.^2 + 2*a2 ;
