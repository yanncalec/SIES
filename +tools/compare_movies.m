function [err, V1i, V2i] = compare_movies(V1, V2)
% err = compare_movies(V1, V2)
% Compare two time series V1 and V2 which are 3D arrays (movie of different sampling rate in
% time) by interpolation. At each instant t, V1(:,:,t) and V2(:,:,t) are 2d images of the same
% dimension, but the length of the two movies are different. We compare
% V1 and V2 by an interpolation in time. The average euclidean distance of each frame is returned.

[m1, n1, l1] = size(V1); 
[m2, n2, l2] = size(V2);

if m1~=m2 || n1~=n2
    error('Two movies must have the same resolution!');
end

if l1 == l2
    err = norm(V1(:)-V2(:), 'fro')/(l1*m1*n1);
    V1i=V1; V2i=V2;
else
    l0 = lcm(l1, l2);
    [xi, yi, zi] = meshgrid(1:m1, 1:n1, 1:l0);

    [x1, y1, z1] = meshgrid(1:m1, 1:n1, linspace(1, l0, l1));
    [x2, y2, z2] = meshgrid(1:m1, 1:n1, linspace(1, l0, l2));

    V1i = interp3(x1,y1,z1,V1,xi,yi,zi);

    V2i = interp3(x2,y2,z2,V2,xi,yi,zi);
    
    err = norm(V1i(:)-V2i(:), 'fro')/(l0*m1*n1);
end
