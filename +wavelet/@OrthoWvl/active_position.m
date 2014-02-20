function [Nx, Ny] = active_position(Z0, Lx, Ly, N1, N2, J, k)
% [Nx, Ny] = active_position(Z0, Lx, Ly, N1, N2, J, k)
%
% Get active positions of the orthogonal wavelet psi or the scaling function phi of a scale. The
% active position is the one where the support of phi_{j,n} or psi_{j,n} intersects the ROI. If
% the qmf index range is in [N1,N2], then the support of 1D phi is [N1,N2] and that of 1D psi is
% [(1-(N2-N1))/2, (1+(N2-N1))/2].
%
% INPUTS:
% Z0: center of the rectangular ROI
% Lx, Ly: width and height of the ROI
% N1, N2: starting and ending index of the qmf, e.g. N1=0, N2=3 for
% Daubechies 2 wavelet.
% J: the scale number
% k: 0,1,2,3 for A-A, A-D, D-A, D-D tensor wavelets. Eg, k=0 corresponds
% to phi, etc.
%
% OUTPUTS:
% Nx: index range in X-axis
% Ny: index range in Y-axis
%
    Ls = N2 - N1; % size of wavelet support

    Ax = [2^(-J)*Z0(1) - 2^(-J-1)*Lx - N2, 2^(-J)*Z0(1) + 2^(-J-1)*Lx - N1];
    Ay = [2^(-J)*Z0(2) - 2^(-J-1)*Ly - N2, 2^(-J)*Z0(2) + 2^(-J-1)*Ly - N1];

    Dx = [2^(-J)*Z0(1) - 2^(-J-1)*Lx - (1+Ls)/2, 2^(-J)*Z0(1) + 2^(-J-1)*Lx - (1-Ls)/2];
    Dy = [2^(-J)*Z0(2) - 2^(-J-1)*Ly - (1+Ls)/2, 2^(-J)*Z0(2) + 2^(-J-1)*Ly - (1-Ls)/2];

    if k==0
        Nx = [ceil(Ax(1)), floor(Ax(2))]; 
        Ny = [ceil(Ay(1)), floor(Ay(2))];
    elseif k==1
        Nx = [ceil(Ax(1)), floor(Ax(2))]; 
        Ny = [ceil(Dy(1)), floor(Dy(2))];
    elseif k==2
        Nx = [ceil(Dx(1)), floor(Dx(2))]; 
        Ny = [ceil(Ay(1)), floor(Ay(2))];
    elseif k==3
        Nx = [ceil(Dx(1)), floor(Dx(2))]; 
        Ny = [ceil(Dy(1)), floor(Dy(2))];
    end
end
