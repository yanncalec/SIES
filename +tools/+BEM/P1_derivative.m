function du = P1_derivative(M,m,L)
% Derivative of the P1 basis (hat function) on a closed curve.
% Inputs:
% M,m: same as in the function P1_basis
% L: the length of the parameterization interval. If the curve is
% parameterized on [0, 2*pi), then L=2*pi

if mod(M, m)
    error('Value Error: the number of boundary points must be a multiplier of the step length.');
end

N = M/m; % the number of hat functions
du = zeros(N,M) ;

dhat = zeros(1,M); 
dhat(1:(2*m+1)) = [ones(1,m), 0, -ones(1,m)];
du(1,:) = circshift(dhat,[0,-m]);

% Then, circular shifting
for i=2:N
    du(i,:) = circshift(du(i-1,:),[0 m])  ;
end

du = M/(L*m) * du';