function [ H, H1 ] = Green2D_Hessian(X, Y)
% Hessian of the 2D Green function G(x) = 1/2/pi * log|x| evaluated at the
% points X-Y. The inputs X, Y are arrays of dimension 2X?, and
% the output H is a matrix whose (m,n) term is the 2X2 Hessian matrix D^2 G(X(:,m) - Y(:,n)).
% and H1=[Dxx Dxy; Dyx Dyy], with Dxx the block matrix of second order
% derivative wrt x taken at X-Y (the tensor matrix), etc.

H = zeros(2*size(X,2), 2*size(Y,2));

% X1 = reshape(X(1,:), [], 1); X2 = reshape(X(2,:), [], 1);
% Y1 = Y(1,:); Y2 = Y(2,:);
% 
% XY1 = repmat(X1, 1, size(Y, 2)) - repmat(Y1, size(X, 2), 1);
% XY2 = repmat(X2, 1, size(Y, 2)) - repmat(Y2, size(X, 2), 1);

XY1 = tools.tensorplus(X(1,:), -Y(1,:));
XY2 = tools.tensorplus(X(2,:), -Y(2,:));

DN = 2*pi* (XY1.^2 + XY2.^2).^2;
M1 = (XY2.^2 - XY1.^2)./DN;
M2 = -2*XY1.*XY2./DN;
M3 = (XY1.^2 - XY2.^2)./DN;

H(1:2:end, 1:2:end) = M1; %Dx^2
H(1:2:end, 2:2:end) = M2; %DxDy
H(2:2:end, 1:2:end) = M2; %DyDx
H(2:2:end, 2:2:end) = M3; %Dy^2

H1 = [[M1 M2]; [M2 M3]];
end

