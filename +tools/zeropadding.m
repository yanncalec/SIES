function [Xz] = zeropadding(X, Nz, mode)
% Zero pad a vector X to a length of Nz.

Nx = length(X);

if nargin<3
    mode = 'tail';
end

if Nz < Nx
    Xz = X;
    %error('Invalid factor n!');
else
    Xz = zeros(Nz, 1); 
    
    if strcmp(mode, 'center')
        %         N0 = floor((Nz-Nx)/2);
        %         Xz(N0+1:N0+Nx) = X;
        if Nx==1
            error('Mode not valid for a scalar');
        else
            if mod(Nx,2) % if Nx odd
                Xz(1:(Nx+1)/2) = X(1:(Nx+1)/2);
                Xz(end-(Nx-1)/2+1:end) = X((Nx+3)/2:end);
            else
                Xz(1:Nx/2+1) = X(1:Nx/2+1); 
                Xz(end-Nx/2+1:end) = X(Nx/2+1:end); % repeat the entry X(N/2+1)
            end
        end
    elseif strcmp(mode, 'head')
        Xz(end-Nx+1:end) = X;
    elseif strcmp(mode, 'tail')
        Xz = zeros(Nz, 1); Xz(1:Nx) = X;
    end
end
