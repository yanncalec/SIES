function [S, R] = make_matrix_SR(cfg, current, Z, impd, order)
% Construct the matrices S and R which are involved in the forward linear operator
% passing from CGPT to SFR data.
% Inputs:
% cfg: acquisition config
% current: surface measurement (the current du/dn)
% Z: reference center
% impd: impedance of the skin
% order: maximum order of the CGPT

    S = zeros(cfg.Ns_total,2*order) ;
    R = cell(1,cfg.Ns_total);

    for s=1:cfg.Ns_total
        Xr = cfg.rcv(s); % s-th receiver
                         
        % The right hand matrix (concerning only the receivers)
        R{s} = -1*PDE.Conductivity_R2.make_matrix_A(Xr, Z, order); 
        
        xs = cfg.src(s); % s-th source
        dipole = cfg.dipole(s); 
        mes = current(s,:); % surface measurement
        
        Omega0 = cfg.Bodies(s); Omega = Omega0.subset(cfg.idxRcv);
        sigma = Omega.sigma;
        normal = Omega.normal;
        
        for m=1:order
            % dipole
            [phim, psim] = phim_psim(m+1, Z(1)-xs(1), Z(2)-xs(2));
            
            A = (-1)^m/2/pi* ( dipole(1)*phim + dipole(2)*psim ) ;
            B = (-1)^m/2/pi* ( dipole(1)*psim - dipole(2)*phim ) ;
            
            % Single layer
            [phim, psim] = phim_psim(m, Xr(1,:) - Z(1), Xr(2,:) - Z(2));
            
            A =  A - 1/2/pi/m * (sigma.*mes) * phim(:) ;
            B =  B - 1/2/pi/m * (sigma.*mes) * psim(:) ;
            
            % Double layer
            [phim, psim] = phim_psim(m+1, Xr(1,:) - Z(1), Xr(2,:) - Z(2));
            
            v1 = phim .* normal(1,:) + psim .* normal(2,:);
            A =  A - impd/2/pi * (sigma.*mes) * v1(:);
            v2 = psim .* normal(1,:) - phim .* normal(2,:);
            B =  B - impd/2/pi * (sigma.*mes) * v2(:);
            
            S(s,2*m-1:2*m) = [A, B];
        end
    end
end

function [phim, psim] = phim_psim(m,x,y)
    [tt,rr] = cart2pol(x, y);
    phim = cos(m*tt) ./ (rr.^m);
    psim = sin(m*tt) ./ (rr.^m);
end
