function dUdn = cylind_wave_Dn(k, m, D)
% Calculate the normal derivative of the cylindrical wave at the boundary D.

    X = D.points;
    Complx = D.cpoints;
    Modul = abs(Complx);
    Angle = angle(Complx);
    
    if 1 %m==1
        if sum(Modul(:)==0) > 0
            error('D contains 0!');
        end
    end
    
    Jn = besselj(m-1, k*Modul);
    J0 = besselj(m,   k*Modul);
    Jp = besselj(m+1, k*Modul);

    t1 = k/2*(Jn - Jp).*X(1,:)./Modul - 1i*m*J0.*X(2,:)./(Modul.^2);
    t2 = k/2*(Jn - Jp).*X(2,:)./Modul + 1i*m*J0.*X(1,:)./(Modul.^2);
    toto = t1.*D.normal(1,:) + t2.*D.normal(2,:);
    dUdn = toto .* exp(1i * m * Angle);
    dUdn = dUdn(:);
end

