function [v, idx] = dist_p2D(z, D)
% [v, idx] = dist_p2D(z, D)
% Euclidean distance between a point z and a domain D.
    
    Dz = D-repmat(z(:), 1, size(D,2));
    [v, idx] = min(sqrt(Dz(1,:).^2 + Dz(2,:).^2));
end