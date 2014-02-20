function Y = func_dh(X, Sm, Rm, M0)
% Y = func_dh(X, Sm, Rm, M0)
% Evaluate the partial derivative of h function of the observation equation in KF
% Inputs:
% X: system state vector
% Sm, Rm: source and receiver matrix returned by matrix_SR_CGPT
% M0: CGPT matrix
% Output:
% Y: partial derivatives
 
% require CGPT M0 be square
if size(M0,1)~=size(M0,2)
    error('CGPT matrix must be square');
end
ord = size(M0,1)/2;

% T0: translation, 2D vector
% Phi0: rotation, scalar
T0 = (X(3) + X(4) * j);
Phi0 = X(5);

% construct the embedding matrix U
u = [1;j];
U = kron(eye(ord), u);
RU = real(U); IU = imag(U);

% Matrix F and its derivatives
% Remark that the partial derivative wrt the position taken at zero is not zero:
% \partial_z F(0) <> 0
% For similar reason, the function dh (and h) cannot be used for tracking only the orientation.

F = zeros(ord); dt1_F = F; dt2_F = F; dp_F = F;
for m=1:ord
    for n=m:ord
        Cmn = nchoosek(n,m); % Lower triangle matrix
        F(m,n) = Cmn * T0^(n-m) * exp(j*m*Phi0);
        if n>m
            % Parital derivative wrt x axis
            dt1_F(m,n) = Cmn * (n-m) * T0^(n-m-1) * exp(j*m*Phi0);
            % Parital derivative wrt y axis
            dt2_F(m,n) = Cmn * (n-m) * T0^(n-m-1) * exp(j*m*Phi0) * j;
        end
        % Parital derivative wrt angle
        dp_F(m,n) = Cmn * T0^(n-m) * exp(j*m*Phi0) * j * m;
    end
end

% Building blocks
RJ = real(U*F); IJ = imag(U*F); % J = U*F;
dt1_RJ = real(U * dt1_F); dt1_IJ = imag(U * dt1_F);
dt2_RJ = real(U * dt2_F); dt2_IJ = imag(U * dt2_F);
dp_RJ  = real(U * dp_F ); dp_IJ  = imag(U * dp_F );

dt1_RMR = dt1_RJ' * M0 * RJ + RJ' * M0 * dt1_RJ;
dt1_RMI = dt1_RJ' * M0 * IJ + RJ' * M0 * dt1_IJ;
dt1_IMR = dt1_IJ' * M0 * RJ + IJ' * M0 * dt1_RJ;
dt1_IMI = dt1_IJ' * M0 * IJ + IJ' * M0 * dt1_IJ;

dt2_RMR = dt2_RJ' * M0 * RJ + RJ' * M0 * dt2_RJ;
dt2_RMI = dt2_RJ' * M0 * IJ + RJ' * M0 * dt2_IJ;
dt2_IMR = dt2_IJ' * M0 * RJ + IJ' * M0 * dt2_RJ;
dt2_IMI = dt2_IJ' * M0 * IJ + IJ' * M0 * dt2_IJ;

dp_RMR = dp_RJ' * M0 * RJ + RJ' * M0 * dp_RJ;
dp_RMI = dp_RJ' * M0 * IJ + RJ' * M0 * dp_IJ;
dp_IMR = dp_IJ' * M0 * RJ + IJ' * M0 * dp_RJ;
dp_IMI = dp_IJ' * M0 * IJ + IJ' * M0 * dp_IJ;

dt1h = RU * dt1_RMR * RU' + RU * dt1_RMI * IU' + IU * dt1_IMR * RU' + IU * dt1_IMI * IU';
dt2h = RU * dt2_RMR * RU' + RU * dt2_RMI * IU' + IU * dt2_IMR * RU' + IU * dt2_IMI * IU';
dph  = RU * dp_RMR * RU' + RU * dp_RMI * IU' + IU * dp_IMR * RU' + IU * dp_IMI * IU';

% Apply then the acquisition operator
Ns=size(Sm,1);
Nr=size(Rm,1);

Y=zeros(Ns*Nr, 5);
Y(:, 3) = reshape(Sm * dt1h * Rm', [], 1);
Y(:, 4) = reshape(Sm * dt2h * Rm', [], 1);
Y(:, 5) = reshape(Sm * dph * Rm', [], 1);
