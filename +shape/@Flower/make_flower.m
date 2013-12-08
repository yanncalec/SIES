function [theta,D,tvec,norm_Tan_square,normal,avec,Sigma,Ksymm] = make_flower(x0,y0,n,k,a,b,phi,epsilon,N)
% This function computes the shape (and its geometrical characteristic) of
% a flower-like set. It is defined as a perturbation of the unit disk.
% Inputs : x0,y0 is the center of the flower
%          a,b are the typical size of the flower
%          phi is an angle parameter
%          n is the number of petals
%          k is the exponent used in the perturbation
%          epsilon is the size of the perturbation
%          N is the number of interpolation points
%          KSymm is the number of rotation periodicity in 2pi

if ~(k>0)
    error 'k must be an stricly positive integer'
end

theta = 2*pi/N*(0:N-1) ;

rot = [[cos(phi) -sin(phi)];[sin(phi) cos(phi)]] ;

% position vector
D = zeros(2,N) ;
% $$$ D(1,:) = x0*ones(1,N) + a*cos(theta).*(ones(1,N)+epsilon*cos(n*theta).^k) ;
% $$$ D(2,:) = y0*ones(1,N) + b*sin(theta).*(ones(1,N)+epsilon*cos(n*theta).^k) ;
% $$$ D = rot*D ;
D(1,:) = cos(theta).*(ones(1,N)+epsilon*cos(n*theta).^k) ;
D(2,:) = sin(theta).*(ones(1,N)+epsilon*cos(n*theta).^k) ;
D = rot*[[a 0];[0 b]]*D + repmat([x0,y0]', 1, N) ;

%velocity vector
tvec = zeros(2,N) ;
tvec(1,:) = -sin(theta).*(ones(1,N)+epsilon*cos(n*theta).^k) - k*epsilon*n*sin(n*theta).*cos(n*theta).^(k-1).*cos(theta) ;
tvec(2,:) = cos(theta).*(ones(1,N)+epsilon*cos(n*theta).^k) - k*epsilon*n*sin(n*theta).*cos(n*theta).^(k-1).*sin(theta) ;
tvec = rot*[[a 0];[0 b]]*tvec ;

norm_Tan_square = tvec(1,:).^2 + tvec(2,:).^2 ;

% unit normal vector
rotation = [[0 1];[-1 0]] ;
normal = rotation*tvec ;
normal = normal./repmat(sqrt(tvec(1,:).^2 + tvec(2,:).^2),2,1) ;


% acceleration vector
avec = zeros(2,N) ;
if(k==1)
    avec(1,:) = -cos(theta).*(ones(1,N)+epsilon*cos(n*theta)) + 2*epsilon*n*sin(theta).*sin(n*theta) - epsilon*n^2*cos(n*theta).* cos(theta);
    avec(2,:) = -sin(theta).*(ones(1,N)+epsilon*cos(n*theta)) - 2*epsilon*n*cos(theta).*sin(n*theta) - epsilon*n^2*cos(n*theta).* sin(theta);
else
avec(1,:) = -cos(theta).*(ones(1,N)+epsilon*cos(n*theta).^k) + 2*k*epsilon*n*sin(theta).*sin(n*theta).*cos(n*theta).^(k-1) + ...
    - k*epsilon*n^2*( cos(n*theta).^k + (k-1)*n*sin(n*theta).^2.*cos(n*theta).^(k-2) ) .* cos(theta);
avec(2,:) = -sin(theta).*(ones(1,N)+epsilon*cos(n*theta).^k) - 2*k*epsilon*n*cos(theta).*sin(n*theta).*cos(n*theta).^(k-1) + ...
    - k*epsilon*n^2*( cos(n*theta).^k + (k-1)*n*sin(n*theta).^2.*cos(n*theta).^(k-2) ) .* sin(theta);
end
avec = rot*[[a 0];[0 b]]*avec ;

% small length element (for integration)
Sigma = 2*pi/N * sqrt(tvec(1,:).^2 + tvec(2,:).^2) ;

% Symmetry of the flower
if a==b
    Ksymm = (mod(k+1,2)+1)*n;
else
    Ksymm = mod(k+1,2)+1;
end



