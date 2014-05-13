function out = make_linop_SCT(cfg, k0, ord)
% Make the linear operator of the acquisition system for the reconstruction of Scattering
% coefficients (SCT).
% INPUTS:
% cfg:
% k0: wave number
% ord: order of SCT
    
    % The matrix related to sources:
    src = cfg.all_src();
    As = PDE.Helmholtz_R2.make_matrix_Src(k0, src, cfg.center, ord);

    % The matrix related to receivers
    if cfg.fixed_rcv
        % If the receiver does not depend on the source, use a simplified model
        % for speed
        Ar = PDE.Helmholtz_R2.make_matrix_Rcv(k0, cfg.rcv(1), cfg.center, ord);
        L = @(x,tflag)tools.linsys.SXR_op(x,As,Ar,tflag); % linear operator
    else
        Ar = {};
        for n=1:cfg.Ns
            Ar{n} = PDE.Helmholtz_R2.make_matrix_Rcv(k0, cfg.rcv(n), cfg.center, ord);
        end
        L = @(x,tflag)tools.linsys.SXR_op_list(x,As,Ar,tflag); % linear operator
    end

    out.L = L;
    out.As = As;
    out.Ar = Ar;

end
