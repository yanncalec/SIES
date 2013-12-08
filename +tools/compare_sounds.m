function [err, V1i, V2i] = compare_sounds(V1, V2)
% err = compare_sounds(V1, V2)
% Compare two time series V1 and V2 which are 1D arrays (sounds of different sampling rate in time)
% by interpolation. We compare V1 and V2 by an interpolation in time. The average euclidean distance
% is returned.

l1 = length(V1); 
l2 = length(V2);

if l1 == l2
    err = norm(V1(:)-V2(:), 'fro')/sqrt(l1);
    % err = norm(V1(:)-V2(:), 'fro');
    V1i=V1; V2i=V2;
else    
    l0 = max(l1, l2);
    % l0 = lcm(l1, l2);
    
    xi = 1:l0;
    x1 = linspace(1, l0, l1);
    x2 = linspace(1, l0, l2);

    V1i = interp1(x1,V1,xi);
    V2i = interp1(x2,V2,xi);
    
    err = norm(V1i(:)-V2i(:), 'fro')/sqrt(l0);
    % err = norm(V1i(:)-V2i(:), 'fro');
end
