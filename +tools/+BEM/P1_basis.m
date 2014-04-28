function u = P1_basis(M,m)
% Compute the set of P1 BEM basis (hat function) on a closed curve whose
% boundary is discretized by M distinct points. A hat function on R is
% defined as: h(t) = 1 - |t| for t in [-1,1]. The discrete hat function is
% defined by 2m+1 points (being 0 at the first and the last points, and 1
% at the middle point). The boundary being periodic, there are M/m basis
% functions, which are equispaced translations of a single hat function.
% The i-th column of the output matrix corresponds to the i-th basis function.
%
% Let X(t) the parameterization of the boundary of D, we define the hat
% function on a boundary as u(X(t)) = h(t). Then:
% \int_{\p D} u(x)v(x) ds(x) = \int_0^{2\pi} h(t) v(X(t)) |X'(t)| dt
% for any function v defined on \p D.
% We will always use the hat function u in the sense above.

if mod(M, m)
    error('Value Error: the number of boundary points must be a multiplier of the step length.');
end

N = M/m; % the number of hat functions
u = zeros(N,M) ;

% We define the hat function centered at the first boundary point, which is
% periodic.
hat = zeros(1,M);
hat(1:(2*m+1)) = [(0:m-1)/m, 1, (m-1:-1:0)/m];
u(1,:) = circshift(hat,[0,-m]);

% Then, circular shifting
for i=2:N
    %     u(i,:) = [u(i-1,M-m:M) u(i-1,1:M-m-1)] ;
    u(i,:) = circshift(u(i-1,:),[0 m])  ;
end

u = u';