function mask = make_WPT_mask(dim, n)
% Make the diagonal mask for selection of AA- or DD-type WPT coefficients.
% INPUTS:
% dim: dimension of the active wavelet coeffs
% n: width of the diagonal mask, n is range of interaction btw wavelets, e.g. n=2
% means each wvl interacts with those inside a square neighborhood of size 5X5.
    
    % x0=tools.extend_matrix(ones(n,n), [N,N]); x0=x0(:)';
    % idx=find(x0);
    % x1 = circshift(x0, [0, -idx(1)]);
    % mask = gallery('circul',x1');
    
    mask = zeros(prod(dim));
    for c=1:dim(2)
        c0=max(1,c-n); c1=min(dim(2),c+n);
        for r=1:dim(1)
            r0=max(1,r-n); r1=min(dim(1),r+n);
            X0=zeros(dim);            
            X0(r0:r1,c0:c1) = 1;
            
            mask((c-1)*dim(1)+r,:) = reshape(X0, [],1);
        end
    end
    
end