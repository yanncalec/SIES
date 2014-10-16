function out = make_linop_CGPT(cfg, ord, symmode)
% Make the linear operator of the acquisition system for reconstruction of CGPT
% INPUTS:
% cfg: configuration of acquisition system
% ord: maximum order of CGPT to be reconstructed
% symmode: true to turn on the symmetry constraint of CGPT matrix

if nargin < 3
	symmode = 0;
end

% The matrix related to sources:
if cfg.nbDirac == 1
	As = PDE.Conductivity_R2.make_matrix_A(cfg.all_src, cfg.center, ord);
else
	src = zeros(2, cfg.nbDirac*cfg.Ns_total);
	
	for s=1:cfg.Ns_total
		src(:, ((s-1)*cfg.nbDirac+1) : s*cfg.nbDirac) = cfg.neutSrc(s); % the positions of diracs of this source satisfying the neutrality condition (see acq.Concerntric.m)
	end
	As0 = PDE.Conductivity_R2.make_matrix_A(src, cfg.center, ord);
	As = kron(eye(cfg.Ns_total),reshape(cfg.neutCoeff, 1, [])) * As0;
end

% The matrix related to receivers
if cfg.Ng == 1
	% If there is only one group, the receiver does not depend on the source and we use a simplified model
	% for speed
	if isa(cfg, 'acq.Coincided') && cfg.nbDirac == 1
		Ar = As;
	else
		Ar = PDE.Conductivity_R2.make_matrix_A(cfg.all_rcv, cfg.center, ord);
	end
	
	if symmode
		L = @(x,tflag)tools.linsys.SXR_op_sym(x,As,Ar,tflag); % linear operator with constraint of symmetry
	else
		L = @(x,tflag)tools.linsys.SXR_op(x,As,Ar,tflag); % linear operator
	end
else
	for n=1:cfg.Ns_total
		Ar{n} = PDE.Conductivity_R2.make_matrix_A(cfg.rcv(n), cfg.center, ord);
	end
	if symmode
		L = @(x,tflag)tools.linsys.SXR_op_list_sym(x,As,Ar,tflag); % linear operator with constraint of symmetry
	else
		L = @(x,tflag)tools.linsys.SXR_op_list(x,As,Ar,tflag); % linear operator
	end
end

out.L = L;
out.As = As;
out.Ar = Ar;
end