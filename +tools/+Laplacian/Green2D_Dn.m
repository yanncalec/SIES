function Gn = Green2D_Dn(X, Y, normal)
% Gn = Green2D_Dn(X, Y, normal)
% Normal derivative of the 2D Green function F(y)=G(y,x) = 1/2/pi * log|y-x| with respect to y on a
% boundary. We compute <DF(y), n_y>. Remark that although mathematically speaking G(x,y)=G(y,x),
% however the order of the arguments must be "X, Y, normal" here.
%
% Inputs:
% X: 2 X M points
% Y, normal: 2 X N boundary points and the normal vector
% Output:
% G is a matrix whose (m, n) term is the normal derivative evlauated at (X(:,m), Y(:,n)).


    [Gx Gy] = tools.Laplacian.Green2D_Grad(Y, X); % Gradient of F(y), which is also the gradient
                                                  % of 1/2/pi*log|x| evaluated at Y-X.
    Gn = diag(normal(1,:))*Gx + diag(normal(2,:))*Gy;
    Gn = Gn';

end
