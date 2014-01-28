function Y = func_h_angle(X, Sm, Rm, M0)
% Y = func_h(X, Sm, Rm, M0)
% Evaluate the h function of the observation equation in KF
% Inputs:
% X: system state vector
% Sm, Rm: source and receiver matrix returned by matrix_SR_CGPT
% M0: CGPT matrix
% Output:
% Y: observation vector (MSR)

% Apply first the translation and rotation on original CGPT
Mc = CGPT_transform(M0, [0,0]', 1, X);

Ns=size(Sm,1);
Nr=size(Rm,1);

Y=zeros(Ns, Nr);

% Apply then the acquisition operator on new CGPT to get the MSR
Y = Sm * Mc * Rm';

% Put the result into vector form
Y=Y(:);