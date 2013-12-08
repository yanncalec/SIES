function [theta,D,tvec,norm_Tan_square,normal,avec,Sigma] = ...
    make_damaged_flower(x0,y0,n,a,b,phi,epsilon,N,tau)
% [theta,D,tvec,norm_Tan_square,normal,avec,Sigma] =  make_damaged_flower(x0,y0,n,a,b,phi,epsilon,N,tau)
% Similar to make_flower, this function gives the shape of a flower
% with the first petal damaged.
%
% Inputs : x0,y0 is the center of the flower
%          a,b are the typical size of the flower
%          phi is an angle parameter
%          n is the number of petals
%          epsilon is the size of the perturbation
%          N is the number of interpolation points
%          tau : tau \in (0,1) is the percentage of damage

theta = 2*pi/N*(0:N-1) ;

% position vector
D = zeros(2,N) ;
% $$$ D(1,:) = x0*ones(1,N) + a*cos(theta).*(ones(1,N)+epsilon*cos(n*theta).^k) ;
% $$$ D(2,:) = y0*ones(1,N) + b*sin(theta).*(ones(1,N)+epsilon*cos(n*theta).^k) ;
% $$$ D = rot*D ;

N0 = ceil(N/n);
theta0 = theta(1:N0); theta1 = theta(N0+1:end);

[~, f, df, ddf] = tools.poly_petal_interp(epsilon, n, 1-tau, theta0);

% figure(); plot(f);
% figure(); plot(df);
% figure(); plot(ddf);

% boundary D
D(1,1:N0) = cos(theta0).*f ;
D(2,1:N0) = sin(theta0).*f ;
D(1,N0+1:end) = cos(theta1).*(ones(1,N-N0)-epsilon*cos(n*theta1)) ;
D(2,N0+1:end) = sin(theta1).*(ones(1,N-N0)-epsilon*cos(n*theta1)) ;

rot = [[cos(phi) -sin(phi)];[sin(phi) cos(phi)]] ;
D = rot*[[a 0];[0 b]]*D + repmat([x0,y0]', 1, N) ;

% velocity vector
tvec = zeros(2,N) ;
tvec(1,1:N0) = -sin(theta0).*f  + cos(theta0).*df;
tvec(2,1:N0) = cos(theta0).*f   + sin(theta0).*df;
tvec(1,N0+1:end) = -sin(theta1).*(ones(1,N-N0)-epsilon*cos(n*theta1)) + epsilon*n*sin(n*theta1).*cos(theta1) ;
tvec(2,N0+1:end) = cos(theta1).*(ones(1,N-N0)-epsilon*cos(n*theta1)) + epsilon*n*sin(n*theta1).*sin(theta1) ;

tvec = rot*[[a 0];[0 b]]*tvec ;

norm_Tan_square = tvec(1,:).^2 + tvec(2,:).^2 ;

% unit normal vector
rotation = [[0 1];[-1 0]] ;
normal = rotation*tvec ;
normal = normal./repmat(sqrt(tvec(1,:).^2 + tvec(2,:).^2),2,1) ;

% acceleration vector
avec = zeros(2,N) ;
avec(1,1:N0) = -cos(theta0).*f - 2*sin(theta0).*df + cos(theta0).*ddf;
avec(2,1:N0) = -sin(theta0).*f + 2*cos(theta0).*df + sin(theta0).*ddf;

avec(1,N0+1:end) =  -cos(theta1).*(ones(1,N-N0)-epsilon*cos(n*theta1)) ...
    - 2*epsilon*n*sin(theta1).*sin(n*theta1) + epsilon*n^2* ...
    cos(n*theta1).* cos(theta1);

avec(2,N0+1:end) = -sin(theta1).*(ones(1,N-N0)-epsilon*cos(n*theta1)) ...
    + 2*epsilon*n*cos(theta1).*sin(n*theta1) + epsilon*n^2* ...
    cos(n*theta1).* sin(theta1);

avec = rot*[[a 0];[0 b]]*avec ;

% small length element (for integration)
Sigma = 2*pi/N * sqrt(tvec(1,:).^2 + tvec(2,:).^2) ;

