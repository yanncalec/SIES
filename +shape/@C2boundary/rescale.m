function [D, tvec, avec, normal] = rescale(D0, theta0, nbPoints, nsize, dspl)
% Compute all variables related to the boundary from the boundary points D0 and the
% parameterization theta0. The new boundary will be reinterpolated with nbPoints, and
% optionally rescaled to fit the rectangle of size
% nsize=[width, height]. The corner singularities can be
% smoothed out by down-sampling (followed by a up-sampling) with the factor dspl (dspl>=1).

% verify the down-sampling factor
if nargin < 5
    dspl = 1;
end

dspl = ceil(dspl);

if dspl<1
    error('Down-sampling factor must be positive!');
end

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

% Re-sampling
nbPoints0 = length(theta0);
idx = [1:dspl:(nbPoints0-1),nbPoints0];
theta = linspace(theta0(1), theta0(end), nbPoints);

% Compute the boundary vectors by cubic-spline interpolation
[D, tvec, avec, normal] = shape.C2boundary.boundary_vec_interpl(D0(:,idx), theta0(idx), theta);
