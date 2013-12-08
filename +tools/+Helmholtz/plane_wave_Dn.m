function [dUdn] = plane_wave_Dn(k, D, Xi)
% [dUdn] = plane_wave_Dn(k, D, Xi)
% Evaluate the normal derivative of the plane wave of given directions Xi at the boundary of D.
% 
% INPUTS:
% k: wave number
% D: an object of C2boundary
% Xi: direction of sources
    
    U = tools.Helmholtz.plane_wave(k, D.points, Xi);

    ddot = tools.tensorprod(D.normal(1,:), Xi(1,:)) + tools.tensorprod(D.normal(2,:), Xi(2,:));
    dUdn = 1i*k*U.*ddot;
end

