function [D, tvec, avec, normal] = rescale_diff(D0, theta0, nbPoints, nsize)
% Compute all variables related to the boundary from the boundary points D0 and the
% parameterization theta0. The new boundary will be reinterpolated with nbPoints, and
% optionally rescaled to fit the rectangle of size
% nsize=[width, height]. The method in use for computing boundary vectors is the finite difference.

% Resize the boundary to the given box
if nargin < 4
    nsize = [];
end

if ~isempty(nsize)
    minx = min(D0(1,:)) ;
    maxx = max(D0(1,:)) ;
    miny = min(D0(2,:)) ;
    maxy = max(D0(2,:)) ;
    
    z0 = [(minx+maxx)/2; (miny+maxy)/2];
    D0 = [(D0(1,:)-z0(1))*(nsize(1)/(maxx-minx)); (D0(2,:)-z0(2))*(nsize(2)/(maxy-miny))];
end

theta = (0:nbPoints-1)/nbPoints*2*pi;
D = interp1(theta0, D0', theta, 'spline');
D = D';
[tvec,avec,normal] = shape.C2boundary.boundary_vec(D, theta);
end
