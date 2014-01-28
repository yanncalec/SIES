function [MSR0, D_L,tvec_L,normal_L,avec_L,Sigma_L] = MSR_simulation_mobile_target(Xs, Xr, kappa, D, tvec, normal, avec, Sigma, Z, Phi)
% [MSR0, D_L,tvec_L,normal_L,avec_L,Sigma_L] = MSR_simulation_mobile_target(Xs, Xr, kappa, D, tvec,
%                                                   ...  normal, avec, Sigma, Z, Phi)
% This function simulates the MSR data of a mobile target, using fixed sources/receivers
% Inputs:
% Xs, Xr: the coordinates of sources and receivers, of dimension 2 X ?
% D: the boundary coordinates of the target, centered at zero, of dimension 2 X ?
% kappa: real conductivity constant
% tvec..Sigma: parameters related to D, see MSR_simulation
% Z: trajectory of mobile target, starting from the origin, of dimension 2 X ?
% Phi: orientation of mobile target, starting from zero
% Outputs:
% MSR0: stream of MSR data corresponding to the trajectory of target, of dimension (Ns X Nr) X length(Z)
% D_L..Sigma_L: a list of transformed shape D and their related boundary parameters
% For recovering the n-th MSR matrix, use: reshape(MSR(:,n), Ns, Nr)

Ns = length(Xs); Nr = length(Xr);
Ntime = length(Z);

MSR0 = zeros(Ns*Nr,Ntime); MSR = MSR0; % MSR data stream;

D_L={}; tvec_L={}; normal_L={}; avec_L={}; Sigma_L={};


theta = 2*pi*(0:length(D)-1)/length(D) ; %boundary parameterization of D

for n=1:Ntime
    [D_L{n},~,tvec_L{n},normal_L{n},avec_L{n},Sigma_L{n}] = boundary_affine_transform(D, theta, tvec, normal, avec, Sigma, Z(:, n), Phi(n), 1);

    toto = MSR_simulation(kappa, D_L{n}, tvec_L{n}, normal_L{n}, avec_L{n}, Sigma_L{n}, Xs, Xr);
    MSR0(:,n) = toto(:);
end

