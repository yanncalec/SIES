function [dU, dG] = make_matrix_gradUG(cfg, Z, fpsi_bg, impd)
% [dU, dG] = make_matrix_gradUG(cfg, current_bg, impd)
% Construct the matrices grad U and grad(dGdn) which are involved in the
% post-processed dipolar expansion eq 4.8
% Inputs:
% cfg: acquisition configuration
% Z: reference point of measurement
% fpsi_bg: the function psi (Be carful, do not use the coefficient vector of
% the P1 basis) solution of the linear system A.2
% impd: impedance
% Outputs:
% dU: a cfg.Ns X 2 matrix. The s-th row correspond to grad U_s(z) with U_s
% the background solution with source at x_s, and z is the center of measurment
% dG: a list of cfg.Nr X 2 matrix, The r-th row correspond to
% grad(dGdn)(x_r, z), with x_r the r-th receiver

    gradSL = zeros(cfg.Ns, 2);
    gradDL = zeros(cfg.Ns, 2);

    [~, H] = tools.Laplacian.Green2D_Hessian(cfg.all_src, Z);
    toto = cfg.all_dipole;
    gradSrc = [diag(toto(1,:)), diag(toto(2,:))] * H; 
    % gradSrc = tools.bdiag(cfg.dipole, 2) * cfg.all_dipole

    for s=1:cfg.Ns
        Fish = cfg.Bodies(s); % Fish's body at s-th position
        
        % gradient of the simple layer potential
        gradSL(s, :) = ops.SingleLayer.eval_grad(Fish, fpsi_bg(:,s), Z);
        
        % gradient of the double layer potential
        gradDL(s, :) = ops.DoubleLayer.eval_grad(Fish, fpsi_bg(:,s), Z);
        
        [~, H] = tools.Laplacian.Green2D_Hessian(cfg.rcv(s), Z);
        DN = [diag(Fish.normal(1, cfg.idxRcv)), diag(Fish.normal(2, cfg.idxRcv))];
        dG{s} = - DN * H;
    end

    dU = gradSrc + gradSL - impd*gradDL;

end
