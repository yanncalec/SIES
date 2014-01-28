function val  = compute_dHdn(obj, sidx)
% Compute the right hand vector dHdn

    if nargin < 2
        sidx = 1:obj.cfg.Ns_total;
    end

    H1 = zeros(obj.nbBEM1, length(sidx)); 
    H2 = zeros(obj.nbIncls*obj.nbBEM2, length(sidx)); 

    for s=1:length(sidx)
        Omega = obj.cfg.Bodies(sidx(s));
        src = obj.cfg.src(sidx(s));
        dipole = obj.cfg.dipole(sidx(s));
        
        % The following need to be checked:
        toto1 = source_vector(Omega, src, dipole);
        H1(:,s) = obj.Psi' * (Omega.sigma(:) .* toto1); % Use P1 basis

        idx = 0;
        for i=1:obj.nbIncls
            % Use P0 basis
            H2(idx+1:idx+obj.nbBEM2, s) = source_vector(obj.D{i}, src, dipole);
            idx = idx+obj.nbBEM2;
        end
    end
    val = [H1; H2];
end

function val  = source_vector(Omega, src, dipole)
% Calculate the scalar product of the source vector dHdn with the
% boundary element basis, involved in eq A.1 and A.5.
% Inputs:
% Omega: the boundary (fish's body or inclusion) on which the source is evaluated
% src: position of the source
% dipole: dipole direction

    [~, hessian] = tools.Laplacian.Green2D_Hessian(Omega.points, src);
    val = [diag(Omega.normal(1,:)) diag(Omega.normal(2,:))] * hessian * dipole;
end
