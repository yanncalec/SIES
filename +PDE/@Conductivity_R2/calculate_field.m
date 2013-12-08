function [F, F_bg, Sx, Sy, mask] = calculate_field(obj, freq, s, z0, width, N)
% Calculate the background potential field and the potential field due to
% the inclusion.
% Inputs:
% freq: the working frequency, a scalar
% s: index of the source
% z0: center of the mesh
% width: width of the mesh
% N: number of points by side
% Outputs:
% F: potential field u
% F_bg: potential field U
% Sx, Sy: coordinates in x and y-axis of the rectangular region
% mask: binary mask where 0 indicates position where the potential may be numerically
% undefined

    epsilon = width/(N-1)/5;

    [Sx, Sy, mask] = obj.D{1}.boundary_off_mask(z0, width, N, epsilon);
    Z = [Sx(:), Sy(:)]';

    for i = 2:obj.nbIncls
        [~, ~, toto] = obj.D{i}.boundary_off_mask(z0, width, N, epsilon);
        mask = mask .* toto; 
    end

    Phi = obj.compute_phi(freq, s);

    % Background field
    Hs = tools.Laplacian.Green2D(Z, obj.cfg.src(s));

    % Compute the perturbation field V=u-Hs by evaluating the single layer potential
    V = reshape(Hs, 1, []);
    for i=1:obj.nbIncls
        V = V + ops.SingleLayer.eval(obj.D{i}, Phi{i}, Z);
    end

    F = reshape(V, [N,N]);
    F_bg = reshape(Hs, [N,N]);

    %% Artificial modification of F is halmful
    % idx = find(1-mask);
    % toto = max(abs(F(find(mask))));
    % F(idx) = toto * sign(F(idx));
end
