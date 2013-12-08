function out = reconstruct_CGPT_analytic(obj, MSR, ord)
% Reconstruct the contracted GPT from data using analytical method.
%
% When the sources and receivers are equally distributed on a circle, use analytical
% formula for reconstruction. The whole linear system is:
% Cs*Ds*CGPT*Dr*Cr' = MSR
% with Cs ~ Ns X K (resp. Cr ~ Nr X K) a matrix depending only on
% source (resp. receivers) angular positions (relative to the center of
% measurement circle), and Ds (resp. Dr) a diagonal matrix depending only
% on Rs (resp. Rr) and K. K is the maximum order of CGPT to
% be recovered. Then we have:
% Cs'*Cs = Ns/2 * Id, if K < Ns/2
% Cr'*Cr = Nr/2 * Id, if K < Nr/2
% The least-square solution is given by:
% X^+ = inv(Ds)*Cs'*MSR*Cr*inv(Dr) /(Ns*Nr/4)

cfg = obj.cfg;
if isa(cfg, 'acq.Concentric') && cfg.equispaced==1 % In this case Ns=Ns_total, Nr=Nr_total
    
    K = min(ord, min(floor((cfg.Ns-1)/2), floor((cfg.Nr-1)/2))) ; % maximum order of CGPT
        
    [Cs, Ds] = make_matrix_CD(cfg.Ns, cfg.radius_src, K);
    [Cr, Dr] = make_matrix_CD(cfg.Nr, cfg.radius_rcv, K);
    out.As = Cs*Ds;
    out.Ar = Cr*Dr;
    
    iDs = diag(1./diag(Ds));
    iDr = diag(1./diag(Dr));
    
    CGPT = 4*iDs*Cs'*MSR*Cr*iDr / cfg.Ns / cfg.Nr;
    out.res = norm(MSR - (out.As * CGPT * out.Ar'), 'fro');
    out.CGPT = CGPT;
    
else
    error('Analytic reconstruction formula can be applied only for equispaced concentric configuration!');        
end
end

function [C,D] = make_matrix_CD(Ns, Rs, order)
% [C,D] = make_matrix_CD(Ns, Rs, order)
% Construct the matrix C and D  (the linear operator L of a equispaced full
% view setting) involved in the model Conductivity_R2 
% Inputs: 
% Ns: number of sources/receivers 
% Rs: radius of the measurement circle
% order: highest order of CGPT

theta = (0:Ns-1)/Ns*2*pi;
tt = (1:order);
T = theta(:) * tt;
C = kron(cos(T), [1,0]) + kron(sin(T), [0,1]);

dd = diag(1./(2*pi*tt.*(Rs.^tt)));
D = kron(dd, eye(2));
end
