function out = reconstruct_SCT_analytic(obj, MSR, freq, ord)

if isa(obj.cfg, 'acq.Planewave') && obj.cfg.equispaced
	
	% ord = min(ord, min(floor((obj.cfg.Ns-1)/2), floor((obj.cfg.Nr-1)/2))) ;
	
	z0 = obj.cfg.center;
	k0 = obj.wavenb_bg(freq);
	
	As = PDE.Helmholtz_R2.make_matrix_Src(k0, obj.cfg.all_src, z0, ord);
	Br = PDE.Helmholtz_R2.make_matrix_Rcv(k0, obj.cfg.rcv(1), z0, ord);
	
	toto = zeros(2*ord+1,1);
	for n=-ord:ord
		toto(n+ord+1) = abs(1/4* besselh(n, 1, k0 * obj.cfg.radius_rcv))^2;
	end
	Dinv = diag(1./toto);
	
	out.SCT = 1/obj.cfg.Ns/obj.cfg.Nr * As'*MSR*Br * Dinv;
	
	out.As = As;
	out.Br = Br;
	
	out.res = norm(MSR - (As*out.SCT*Br'), 'fro');
else
	error(['Analytic reconstruction formula can be applied only for equispaced circular planewave ' ...
		'configuration!']);
end

