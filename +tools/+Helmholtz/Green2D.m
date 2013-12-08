function [Ga] = Green2D(k, X, Y)
% 2D Green function for Helmholtz equation (\Delta+k^2) f = 0.
%
% \Gamma^k(x) = -sqrt(-1)/4*besselh(0,k*|X-Y|) evaluated at the points X-Y, and k is the wave
% number. The inputs X, Y are arrays of dimension 2X?, and the output Ga is a matrix whose (m,n)
% term is \Gamma^k evluated at (X(:,m) - Y(:,n)).

    if size(X,1)~=2 || size(Y,1)~=2
        error('The inputs X and Y must have two rows!');
    end
    % X = reshape(X, 2, []); 
    % Y = reshape(Y, 2, []);

    XY1 = tools.tensorplus(X(1,:), -Y(1,:));
    XY2 = tools.tensorplus(X(2,:), -Y(2,:));

    SDN = sqrt(XY1.^2 + XY2.^2);
    Ga = -1i/4*besselh(0,1,k*SDN);

end
