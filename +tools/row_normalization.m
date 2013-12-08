function [Asn, toto] = row_normalization(As)
% Normalize the row vectors of a matrix such that each row has euclidian norm 1

    toto = sqrt(sum(As.*As, 2));
    idx = find(toto>0); 
    dd = toto; 
    dd(idx) = 1./toto(idx);
    Asn = diag(dd) * As; % normalize each row
end
