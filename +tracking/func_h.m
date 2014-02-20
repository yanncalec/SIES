function Y = func_h(X, Sm, Rm, M0)
% Y = func_h(X, Sm, Rm, M0)
% Evaluate the h function of the observation equation in EKF
% Inputs:
% X: system state vector
% Sm, Rm: measurement matrices correspondind to sources and receivers
% M0: CGPT matrix
% Output:
% Y: observation vector (MSR)

T0 = X(3:4); 
Phi0 = X(5);

% Apply first the translation and rotation on original CGPT
Mc = asymp.CGPT.CGPT_transform(M0, T0, 1, Phi0);

% Apply then the acquisition operator on new CGPT to get the MSR
Y = Sm * Mc * Rm';

% Put the result into vector form
Y=Y(:);
