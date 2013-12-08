function W = theoretical_SCT(D, ord, freq, pmtt_bg, pmeb_bg)
% Theoretical value of scattering coefficients

    k0 = freq*sqrt(pmeb_bg*pmtt_bg);

    P = zeros(2*ord+1, D.nbPoints);
    C = zeros(D.nbPoints, 2*ord+1);

    matrix_A = PDE.Helmholtz_R2.system_matrix_block_A11(freq, D, 'P0', 1);
    matrix_B = PDE.Helmholtz_R2.system_matrix_block_A12(freq, D, 'P0', 1, pmeb_bg, pmtt_bg);
    matrix_C = PDE.Helmholtz_R2.system_matrix_block_A21(freq, D, 'P0', 1);
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

    %    W = W.';
end