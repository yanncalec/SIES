function W = theoretical_SCT(D, pmeb, pmtt, pmeb_bg, pmtt_bg, ord, freq)
% Calculate the scattering coefficients of (multiple) shape(s) D:
%
%  W_{mn}[D] = \sum_{l=1}^L \int_{\p D_l} C_n(y)^* \psi_m(y) ds(y)
%
% with C_n(y) the cylindrical wave and \psi_m the solution to the
% following system:
%
% S_{D}^{kl}[\phi_l] - S_{D}^{k0}[\psi] = C_m   on \p D,
% \mu_*^{-1} ddn S_{D}^{k}[\phi]|_ - \mu_0^{-1} ddn
% S_{D}^{k0}[\psi]|_+ = \mu_0^{-1} ddn C_m
%
% Inputs:
% D: a C2boundary object
% pmtt, pmeb: permittivity and permeability of inclusion
% pmtt_bg, pmeb_bg: permittivity and permeability of the background
% ord: maximum order of SCT to be computed
% freq: list of frequency

if length(D)>1
    erros('Multi-inclusions are not suppoerted in the current version!');
end

if iscell(D)
    D = D{1};
end

k0 = freq*sqrt(pmeb_bg*pmtt_bg);

P = zeros(2*ord+1, D.nbPoints);
C = zeros(D.nbPoints, 2*ord+1);

matrix_A = PDE.Helmholtz_R2.system_matrix_block_A11(freq, D, 'P0', 1, pmeb, pmtt);
matrix_B = PDE.Helmholtz_R2.system_matrix_block_A12(freq, D, 'P0', 1, pmeb_bg, pmtt_bg);
matrix_C = PDE.Helmholtz_R2.system_matrix_block_A21(freq, D, 'P0', 1, pmeb, pmtt);
matrix_D = PDE.Helmholtz_R2.system_matrix_block_A22(freq, D, 'P0', 1, pmeb_bg, pmtt_bg);
matrix_BEM = [[matrix_A matrix_B] ; [matrix_C matrix_D]] ;

for m=-ord:ord
    b = PDE.Helmholtz_R2.source_vector_cw(D, freq, m, pmeb_bg, pmtt_bg);
    sol = matrix_BEM\b;
    vpsi = sol(size(sol,1)/2+1:end,:);
    
    P(m+ord+1, :) = reshape(vpsi, 1, []);
end

for n=-ord:ord
    C(:, n+ord+1) = (-1)^n * tools.Helmholtz.cylind_wave(k0, -n, D.points);
end

W = P*diag(D.sigma)*C; % definition in 5.15

end
