function [F, F_bg, Sx, Sy, mask] = calculate_field(obj, f, s, z0, width, N, fpsi_bg, fpsi, fphi)
% Calculate the background potential field and the potential field due to
% the inclusion
% Inputs:
% vpsi_bg: the coefficient vector of background eq (A.2)
% vpsi, vphi: the coefficient vector of eq (A.5)
% f: index of the frequency
% s: index of the source
% xlim: interval in x-axis of rectangular region
% ylim: interval in y-axis of rectangular region
% dh: sampling step
% Outputs:
% F: potential field u
% F_bg: potential field U
% SX, SY: coordinates in x and y-axis of the rectangular region

    Omega = obj.cfg.Bodies(s);
    src = obj.cfg.src(s);
    dipole = obj.cfg.dipole(s);
    
    epsilon = width/(N-1)/5; % user specified constant for the mask
    
    [Sx, Sy, mask] = Omega.boundary_off_mask(z0, width, N, epsilon);
    Z = [Sx(:), Sy(:)]';

    for i = 1:obj.nbIncls
        [~, ~, toto] = obj.D{i}.boundary_off_mask(z0, width, N, epsilon);
        mask = mask .* toto; 
    end
    
    [Gx, Gy] = tools.Laplacian.Green2D_Grad(Z, src);
    Hs = reshape(Gx*dipole(1)+Gy*dipole(2), 1, []);
    
    fpsi_bg = fpsi_bg(:, s);
    fpsi = fpsi{f}(:, s);
    fphi = squeeze(fphi{f}(:, s, :));
    
    %% For total field
    % Evaluate the simple and double layer potentials
    Ss = ops.SingleLayer.eval(Omega, fpsi, Z);
    Ds = ops.DoubleLayer.eval(Omega, fpsi, Z);
    
    V = Hs + Ss - obj.impd * Ds ;

    for i=1:obj.nbIncls
        % contribution of the i^th inclusion 
        V = V + ops.SingleLayer.eval(obj.D{i}, fphi(:,i), Z);
    end
    F = reshape(V, [N,N]);
    
    %% For background field
    % Evaluate the simple and double layer potentials
    Ss = ops.SingleLayer.eval(Omega, fpsi_bg, Z);
    Ds = ops.DoubleLayer.eval(Omega, fpsi_bg, Z);
    
    V = Hs + Ss - obj.impd * Ds;
    F_bg = reshape(V, [N,N]);
end
