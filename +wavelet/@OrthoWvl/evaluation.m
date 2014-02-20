function [val, dX, dY] = evaluation(X, J, N, Tx, Ty, Gx, Gy)
% [val, dX, dY] = evaluation(X, J, N, Tx, Ty, Gx, Gy)
%
% Evaluate by looking-up table, the tensor-product scaling or wavelet functions of scale J and
% indexes N at the coordinates X, as well as their partial derivatives. More precisely, the function
% to be evaluated is 2^(-J) psi(2^(-J)(X-2^J*N)), with J an integer, X and N are 2-by-? arrays and N
% takes integer values. The 1D tables Tx and Ty contain precalculated values of the 1D mother
% function psi_1(x) and psi_2(y) on the uniform 1D lattice Gx and Gy
% respectively. psi(x,y)=psi_1(x)*psi_2(y) is compactly supported and its support is contained in Gx
% X Gy.
%
% OUTPUT:
% val: a |X|-by-|N| matrix with each column the values of psi_{j,n}
% evaluated at X.
% dX, dY: similar to val but for partial derivatives.

    dx = Gx(2) - Gx(1); % sampling step of the grid Gx
    dy = Gy(2) - Gy(1); % sampling step of the grid Gy
                     
    % dh = (G(end)-G(1)) / (length(G)-1); 

    % Make the derivative table    
    % dTx = diff(Tx) / dx; dTy = diff(Ty) / dy;  
    % dTx = [dTx 0]; dTy = [dTy 0];

    dTx = tools.mdiff(Tx)/dx; dTy = tools.mdiff(Ty)/dy;

    % Lf=100;
    % dTx = cconv(dTx,ones(1,Lf), length(Tx))/Lf;
    % dTy = cconv(dTy,ones(1,Lf), length(Ty))/Lf;
    
    %% Looking-up table 

    % C version
    [Vx, dVx] = wavelet.OrthoWvl.lookup_table_wvl_mex(2^(-J)*X(1,:), -N(1,:), Gx(1), Gx(end), dx, Tx, dTx);
    [Vy, dVy] = wavelet.OrthoWvl.lookup_table_wvl_mex(2^(-J)*X(2,:), -N(2,:), Gy(1), Gy(end), dy, Ty, dTy);

    % % Matlab version:
    % Cx = tools.tensorplus(2^(-J)*X(1,:), -N(1,:)); % True coordinates of evaluation
    % Cy = tools.tensorplus(2^(-J)*X(2,:), -N(2,:));

    % Mx = (Cx>=Gx(1)) .* (Cx<=Gx(end)); % mask
    % Ix = find(Mx); % index
    % Kx = floor((Cx(Ix)-Gx(1))/dx) + 1; % conversion from coordinate to index
    % Vx = zeros(size(Cx)); Vx(Ix) = Tx(Kx);
    % dVx = zeros(size(Cx)); dVx(Ix) = dTx(Kx);

    % My = (Cy>=Gy(1)) .* (Cy<=Gy(end)); % mask
    % Iy = find(My); % index
    % Ky = floor((Cy(Iy)-Gy(1))/dy) + 1; % conversion from coordinate to index
    % Vy = zeros(size(Cy)); Vy(Iy) = Ty(Ky);
    % dVy = zeros(size(Cy)); dVy(Iy) = dTy(Ky);

    % clear Mx Ix Kx My Iy Ky;

    %% From 1d value to 2d value by tensor product
    val = 2^(-J) * (Vx.*Vy); 
    dX = 2^(-2*J) * (dVx.*Vy);
    dY = 2^(-2*J) * (Vx.*dVy);

end

