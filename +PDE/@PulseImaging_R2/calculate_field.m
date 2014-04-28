function [F, F_bg, Sx, Sy, mask] = calculate_field(obj, Ntime, s, z0, width, N)
% Calculate the background potential field and the potential field due to
% the inclusion.
% Inputs:
% Tmax: time interval [0, Tmax]
% s: index of the source
% z0: center of the mesh
% width: width of the mesh
% N: number of points by side
% Outputs:
% F: potential field u (time-dependent)
% F_bg: potential field U  (time-dependent)
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

    Phi = obj.compute_phi(Ntime, s);

    Hs = reshape(obj.cfg.neutCoeff, 1, []) * tools.Laplacian.Green2D(obj.cfg.neutSrc(s), Z);

    F = {}; F_bg = {};
    for t=1:Ntime
        V = Hs * obj.waveform(t);
        % Background field
        F_bg{t} = reshape(V, [N,N]);

        for i=1:obj.nbIncls
            V = V + ops.SingleLayer.eval(obj.D{i}, squeeze(Phi{t}(:, i, :)), Z);
        end
        % total field
        F{t} = reshape(V, [N,N]);
    end
    
end
