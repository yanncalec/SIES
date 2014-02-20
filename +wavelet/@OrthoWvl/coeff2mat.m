function M = coeff2mat(obj, C)
% Put the coefficient vector into the classical wavelet tree matrix form

    [A,D] = obj.vec2scl(C);
    nlevel = size(D,1);
    M = zeros(size(obj.wmask.A) * 2^nlevel);
    
    rowdim = size(obj.wmask.A, 1); coldim = size(obj.wmask.A, 2); % base dimension
    M(1:rowdim, 1:coldim) = obj.extract_coeffs(C, 0);

    rp0 = rowdim+1; cp0 = coldim+1; % row and column pointer

    for j=1:nlevel
        rp = rp0; cp = cp0-coldim; 
        C0 = obj.extract_coeffs(C, j, 1); % norm(C0,'fro')
        toto = tools.extend_matrix(C0, [rowdim, coldim]); % norm(toto,'fro')
        M(rp:rp+rowdim-1, cp:cp+coldim-1) = toto;

        rp = rp0-rowdim; cp = cp0;
        C0 = obj.extract_coeffs(C, j, 2); % norm(C0,'fro')
        toto = tools.extend_matrix(C0, [rowdim, coldim]); % norm(toto,'fro')
        M(rp:rp+rowdim-1, cp:cp+coldim-1) = toto;
        
        rp = rp0; cp = cp0;
        C0 = obj.extract_coeffs(C, j, 3); % norm(C0,'fro')
        toto = tools.extend_matrix(C0, [rowdim, coldim]); %norm(toto,'fro')
        M(rp:rp+rowdim-1, cp:cp+coldim-1) = toto;

        rp0 = rp0 + rowdim; cp0 = cp0 + coldim;
        rowdim = rowdim*2; coldim = coldim*2;

        % M(rp:rp+rowdim-1, cp:cp+coldim-1) = wextend(2, 'zpd', C0, [rowdim ...
        %                     - cdim(j,1), coldim - cdim(j,2)]);
    end
end
