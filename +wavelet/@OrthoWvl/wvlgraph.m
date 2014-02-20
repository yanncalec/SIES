function [Tphi, Tpsi, Gphi, Gpsi, Tdphi, Tdpsi] = wvlgraph(hcoeff, gcoeff, niter)
% [Tphi, Tpsi, Gphi, Gpsi] = wvlgraph(hcoeff, niter)
% Calculate the graph of 1D orthogonal wavelets using an iterative procedure of convolution and up-sampling.
% INPUTS:
% order: 2...10, the order of wavelet
% niter: number of iteration
% OUTPUTS:
% Tphi: table of values phi(x) of the scaling function phi on its support
% Tpsi: table of values psi(x) of the mother wavelet psi on its support
% Gphi: the support of phi
% Gpsi: the support of psi
    
    if nargin < 2
        niter = 10;
    end

    L = length(hcoeff)-1; % support of the scaling function is in [0, L]    
    fct = niter/2;
    %% Calculate the scaling function phi
    C0=hcoeff;
    for n=1:niter-1
        C = zeros(2*length(C0)-1,1);
        C(1:2:end)=C0;
        C0 = conv(C,hcoeff,'full');
    end
    Tphi = reshape(C0 * 2^fct, 1, []);
    
    %% Calculate the mother wavelet psi
    C0 = gcoeff;
    for n=1:niter-1
        C = zeros(2*length(C0)-1,1);
        C(1:2:end)=C0;
        C0 = conv(C,hcoeff,'full');
    end
    Tpsi = reshape(C0 * 2^fct, 1, []);
    
    %% Calculate the derivative of scaling function phi
    C0=hcoeff;
    for n=1:niter-1
        C = zeros(2*length(C0)-1,1);
        C(1:2:end)=C0;
        C0 = conv(C,2*hcoeff,'full');
    end
    Tdphi = reshape(C0 * 2^-fct, 1, []);
    
    %% Calculate the derivative of mother wavelet psi
    C0=gcoeff;
    for n=1:niter-1
        C = zeros(2*length(C0)-1,1);
        C(1:2:end)=C0;
        C0 = conv(C,2*hcoeff,'full');
    end
    Tdpsi = reshape(C0 * 2^-fct, 1, []);

    Gphi = linspace(0, L, length(Tphi)); 
    Gpsi = linspace((1-L)/2, (1+L)/2, length(Tpsi));
end
