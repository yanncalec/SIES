function [Asn, toto] = col_normalization(As)
% Normalize the column vectors of a matrix such that each column has euclidian norm 1

    toto = sqrt(sum(As.*As, 1));
    idx = find(toto>0); 
    dd = toto; 
    dd(idx) = 1./toto(idx);
    Asn = As * diag(dd); % normalize each column
end
