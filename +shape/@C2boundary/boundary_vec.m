function [tvec,avec,normal] = boundary_vec(D, t)
% Given a curve D (the boundary of a simply connected region) and its
% parameterization t, evaluate the following values:
% tvec: tangent vector
% normal: outward normal vector
% avec: acceleration vector
%
% Remark that the first and the last elements in D and t must NOT
% be the same (not tired-off)

% Verify that D has dimension 2 X N

% tangent vector
tvec=curve_derivative(D,t);

% outward normal vector
rotation = [[0 1];[-1 0]] ;
normal = rotation*tvec ;
normal = normal./repmat(sqrt(normal(1,:).^2+normal(2,:).^2),2,1) ;

% acceleration vector
avec = curve_derivative2(D, t);
%Sigma = 2*pi/length(t) * sqrt(tvec(1,:).^2 + tvec(2,:).^2) ; % length of small line element ( for integration )

function tvec = curve_derivative(D, t)
% tvec = curve_derivative(D, t)
%
% First order derivative of a curve (tangent vector).
% D: boundary point coordinates. 2 X N array, D(:,n) is a 2d vector.
% t: parameterization of D. t takes its value in [0, 2*pi).
%
% D is parameterized as [r(t)*cos(t), r(t)*sin(t)]. It's a circuit but
% with D(:,1)~=D(:,end).

dD = (D - circshift(D,[0,1]));
dt = mod(t-circshift(t,[0,1]), 2*pi);

tvec(1,:) = dD(1,:)./dt;
tvec(2,:) = dD(2,:)./dt;

function avec = curve_derivative2(D, t)
% avec = curve_derivative2(D, t)
% Second order derivative of a closed curve (curvature vector).
%
% Remark that the first and the last elements in D and t must NOT
% be the same (no superposition)

Dr=circshift(D,[0,1]);
Dl=circshift(D,[0,-1]);
dD = Dr-2*D+Dl;

%dt=(2*pi/max(size(D))^2;

dt2 = mod(t-circshift(t,[0,1]), 2*pi).^2;
avec(1,:) = dD(1,:)./dt2;
avec(2,:) = dD(2,:)./dt2;
avec=circshift(avec,[0,-1]);            % correct the offset
