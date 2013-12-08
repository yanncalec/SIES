function [Gx Gy] = Green2D_Grad(X, Y)
% Gradient of the 2D Green function G(x) = 1/2/pi * log|x| evaluated at the
% points X-Y. The inputs X, Y are arrays of dimension 2X?, and the output
% Gx (resp. Gy) is a matrix whose (m,n) term is the derivative Dx (resp.
% Dy) of G evlauated at (X(:,m) - Y(:,n)).

    XY1 = tools.tensorplus(X(1,:), -Y(1,:));
    XY2 = tools.tensorplus(X(2,:), -Y(2,:));

    DN = XY1.^2 + XY2.^2;
    Gx = 1/2/pi * XY1 ./ DN;
    Gy = 1/2/pi * XY2 ./ DN;

end
