function [ M ] = rectangular_mesh(Z, width, height, N)
% Generation of a mesh on a rectangular region

sx = linspace(Z(1)-width/2, Z(1)+width/2, N);
sy = linspace(Z(2)-height/2, Z(2)+height/2, N);
[SX, SY] = meshgrid(sx, sy); 
M=[SX(:) SY(:)]';

end

