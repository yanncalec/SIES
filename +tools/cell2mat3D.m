function Xm = cell2mat3D(Xc)
% Convert a cell of 2D matrix to a 3D matrix:
% Inputs:
% Xc: a cell, Xc{t} is a matrix of dim N1 X N2
% Output:
% Xm: a matrix of dim N1 X N2 X T, with Xm(:,:,t) = Xc{t}

T = length(Xc);

[N1,N2] = size(Xc{1});

Xm = zeros(N1,N2,T);

for t=1:T
    Xm(:,:,t) = Xc{t};
end
