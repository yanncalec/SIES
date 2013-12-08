function [Gx Gy] = Green2D_Grad(k, X, Y)
% Gradient of the 2D Green function of the Helmholtz equation evaluated at the points X-Y.
% The Green function \Gamma^k(x)=-i/4*H_0^1(k|x|)
% The derivative of the hankel function H_0^1(t) = J_0(t)+i*Y_0(t) is -H_1^1(t), because
% J_0'(t)=-J_1(t), Y_0'(t)=-Y_1(t)
    
    XY1 = tools.tensorplus(X(1,:), -Y(1,:));
    XY2 = tools.tensorplus(X(2,:), -Y(2,:));

    SDN = sqrt(XY1.^2 + XY2.^2);

    toto = besselh(1,1,k*SDN);
    Gx = 1i*k/4*toto*XY1./SDN;
    Gy = 1i*k/4*toto*XY2./SDN;
end
