function V = matrix_norm_by_order(X, omode)

if nargin < 2
    omode = 'normal';
end

[M,N] = size(X);
if M~=N
    error('Matrix must be square!');
end

if strcmp(omode, 'center')
    ord = floor(N/2);
    V = zeros(1, ord);
    
    if 2*ord == N
        rr = ord; 
        cc = ord;

        for K=0:ord-1
            toto = X((rr-K):(rr+K+1), (cc-K):(cc+K+1));
            V(K+1) = norm(toto, 'fro');
        end
    else
        rr = ord + 1;
        cc = ord + 1;        

        for K=1:ord
            toto = X((rr-K):(rr+K), (cc-K):(cc+K));
            V(K) = norm(toto, 'fro');
        end
        V = [norm(X(rr,cc), 'fro') V];
    end        
else
    V = zeros(1, N);
    for K=1:N
        V(K) = norm(X(1:K, 1:K), 'fro');
    end
end
