function W = wvl_vec2scl(C, L0)

    P=cumsum(L0(1:end-1))+1;
    W=[];
    W{1} = C(1:P(1)-1);

    for n = 1:length(P)-1
        W{n+1} = C(P(n):P(n+1)-1);
    end

    % norm(cell2mat(W)-X)
end