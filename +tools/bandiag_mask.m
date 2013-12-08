function M = bandiag_mask(N, bw)
% M = bandiag_mask(N, bw)
% Return a band diagonal mask of size N-by-N of the width bw. bw must be odd.

    if N<bw
        error('Dimension error!');
    end
    
    aa=zeros(1, N);
    w = max(1, floor((bw+1)/2));
    aa(1:w)=1;
    M = toeplitz(aa, aa');
end