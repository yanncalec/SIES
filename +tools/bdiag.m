function [ D ] = bdiag(A, idx)
% [ D ] = bdiag(A, idx)
% Make a blockwise diagonal matrix using the row or the column of a matix.
% INPUTS:
% A: input matrix
% idx: 1 for row-wise diagonal matrix, otherwise for column-wise diagonal matrix
[M,N] = size(A);

if idx==1
    D = zeros(M*N, N);
    for n = 1:N
        D(((n-1)*M+1):(n*M), n) = A(:,n);        
    end
else
    D = zeros(M, M*N);
    for m = 1:M
        D(m, ((m-1)*N+1):(m*N)) = A(m,:);        
    end    
end

end

