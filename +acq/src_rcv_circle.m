function [Xs, Thetas, Xscell] = src_rcv_circle(Na, N0, R0, Z, theta, aov)
% Generate sources/receivers placed on concentric arcs of radius
% R0 centered at a point Z. There are in total Na arcs which are equally
% distributed between the angular range [0, aov), and each has the same angular
% coverage which is theta.
%
% Inputs:
% -Na: number of arcs
% -N0: number of sources/receivers per arc
% -R0: radius of measurement circle
% -Z: center of measurement circle
% -theta: angle of each arc
% -aov: (optional) total angle of view, 2*pi by default
%
% Outputs:
% -Xs: the positions of sources/receivers, dimension 2 X (N0*Na)
% -Thetas: angle of each source/receiver
% -Xscell: Xs in cell format. Xscell{n} is the Xs of the n-th arc

if nargin<6
	aov=2*pi;
end

Xs=zeros(2,N0*Na);
Thetas=zeros(1,N0*Na);
%theta = min(theta, aov/Na);

for n=1:Na
	tt0 = (n-1)/Na*aov;
	tt = tt0 + (0:N0-1)/N0 * theta;
	
	rr = R0 * [cos(tt); sin(tt)];
	Xs(:, ((n-1)*N0+1):(n*N0)) = rr;
	Thetas(((n-1)*N0+1):(n*N0)) = tt;
	
	Xscell{n} = Xs(:, ((n-1)*N0+1):(n*N0)) + repmat(Z(:), 1, size(rr,2));
end

Xs = Xs + repmat(Z(:), 1, size(Xs,2));

end