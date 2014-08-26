function [points,tvec,avec,normal] = boundary_vec_interpl(points0, theta0, theta)
% Given a curve (the boundary of a simply connected region) and its
% parameterization, evaluate the tangent, acceleration and normal vector 
% by downsampling then resampling (interpolation by cubic-spline)
%
% Inputs:
% points0: 2 X N array, the coordinates of the curve
% theta: parameterization, eg, 0..N-1
% dspl: downsampling factor
% 
% Outputs:
% points: resampled curve 
% tvec: tangent vector
% normal: outward normal vector
% avec: acceleration vector
%
% This function is useful when dealing with the singularities, and it is
% more accurate than boundary_vec.m.
%
% Remark that the first and the last elements in D and t must NOT
% be the same (not tired-off)

if nargin < 3
    theta = theta0;
end

fx = csapi(theta0(:), points0(1,:));
% fx = csape(theta0(:)', points0(1,:), 'periodic');
px = fnval(fx, theta);

fy = csapi(theta0(:), points0(2,:));
% fy = csape(theta0(:)', points0(2,:), 'periodic');
py = fnval(fy, theta);

points = [reshape(px,1,[]); reshape(py,1,[])];

dfx = fnder(fx,1);
dfy = fnder(fy,1);
tx = fnval(dfx, theta);
ty = fnval(dfy, theta);
tvec = [reshape(tx,1,[]); reshape(ty,1,[])];

rotation = [[0 1];[-1 0]] ;
normal = rotation*tvec ;

normal = normal./repmat(sqrt(normal(1,:).^2+normal(2,:).^2),2,1);

ddfx = fnder(fx,2);
ddfy = fnder(fy,2);
ax = fnval(ddfx, theta);
ay = fnval(ddfy, theta);
avec = [reshape(ax,1,[]); reshape(ay,1,[])];
