function out = reconstruct_CGPT_analytic(obj, MSR, ord)
% Reconstruct the contracted GPT from data using analytical method.
%
% The analytical method of reconstruction is applicable only when the sources and
% receivers are both equally distributed on a circle. The whole linear system is:
%
% Cs*Ds*CGPT*Dr*Cr' = MSR
%
% with Cs ~ Ns X K (resp. Cr ~ Nr X K) a matrix depending only on
% source (resp. receivers) angular positions (relative to the center of
% measurement circle), and Ds (resp. Dr) a diagonal matrix depending only
% on Rs (resp. Rr) and K. K is the maximum order of CGPT to
% be recovered. Then we have:
% Cs'*Cs = Ns/2 * Id, if K < Ns/2
% Cr'*Cr = Nr/2 * Id, if K < Nr/2
% The least-square solution is given by:
%
% X^+ = inv(Ds)*Cs'*MSR*Cr*inv(Dr) /(Ns*Nr/4)
%
% Inputs:
% MSR: the MSR data matrix, if it is time-dependent or multifrequency, MSR
% need to be a cell.
% ord: maximum order of reconstruction

if ~iscell(MSR) % Convert to a cell, for compability with PulseImaging_R2 class
	MSR = {MSR};
end

cfg = obj.cfg;
if isa(cfg, 'acq.Concentric') && cfg.equispaced==1 % In this case Ns=Ns_total, Nr=Nr_total
	if cfg.nbDirac == 1
		K = min(ord, min(floor((cfg.Ns-1)/2), floor((cfg.Nr-1)/2))) ; % maximum order of CGPT
		% K = ord;
		
		[Cs, Ds] = make_matrix_CD(cfg.Ns, cfg.radius_src, K);
		[Cr, Dr] = make_matrix_CD(cfg.Nr, cfg.radius_rcv, K);
		
		iDs = diag(1./diag(Ds));
		iDr = diag(1./diag(Dr));
		out.As = Cs*Ds; out.Ar = Cr*Dr;
		
		for t=1:length(MSR);
			out.CGPT{t} = 4*iDs*Cs'*MSR{t}*Cr*iDr / cfg.Ns / cfg.Nr;
			out.res{t} = norm(MSR{t} - (out.As * out.CGPT{t} * out.Ar'), 'fro');
		end
	else
		src = zeros(2, cfg.nbDirac*cfg.Ns_total);
		
		for s=1:cfg.Ns_total
			src(:, ((s-1)*cfg.nbDirac+1) : s*cfg.nbDirac) = cfg.neutSrc(s); % the positions of diracs of this source satisfying the neutrality condition (see acq.Concerntric.m)
		end
		As0 = PDE.Conductivity_R2.make_matrix_A(src, cfg.center, ord);
		out.As = kron(eye(cfg.Ns_total),reshape(cfg.neutCoeff, 1, [])) * As0;
		out.Ar = PDE.Conductivity_R2.make_matrix_A(cfg.all_rcv, cfg.center, ord);
		
		for t=1:length(MSR);
			out.CGPT{t} = pinv(out.As) * MSR{t} * pinv(out.Ar');
			out.res{t} = norm(MSR{t} - (out.As * out.CGPT{t} * out.Ar'), 'fro');
		end
	end
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
