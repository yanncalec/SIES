function Y = func_dh_angle(X, Sm, Rm, M0)
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

Phi0 = X;

% construct the embedding matrix U
u = [1;j];
U = kron(eye(ord), u);
RU = real(U); IU = imag(U);

F = zeros(ord); dp_F = F;
for m=1:ord
    F(m,m) = exp(j*m*Phi0);
    dp_F(m,m) = exp(j*m*Phi0) * j * m;
end

RJ = real(U*F); IJ = imag(U*F); % J = U*F;
dp_RJ  = real(U * dp_F ); dp_IJ  = imag(U * dp_F );

dp_RMR = dp_RJ' * M0 * RJ + RJ' * M0 * dp_RJ;
dp_RMI = dp_RJ' * M0 * IJ + RJ' * M0 * dp_IJ;
dp_IMR = dp_IJ' * M0 * RJ + IJ' * M0 * dp_RJ;
dp_IMI = dp_IJ' * M0 * IJ + IJ' * M0 * dp_IJ;

dph  = RU * dp_RMR * RU' + RU * dp_RMI * IU' + IU * dp_IMR * RU' + IU * dp_IMI * IU';

% Apply then the acquisition operator
Ns=size(Sm,1);
Nr=size(Rm,1);

Y=zeros(Ns*Nr, 1);
Y = reshape(Sm * dph * Rm', [], 1);
