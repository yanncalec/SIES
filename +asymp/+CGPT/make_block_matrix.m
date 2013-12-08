function KsdS = make_block_matrix(D)
% Construct some blocks of the matrix A in the system A[phi]=b, where phi is the auxiliary
% function(s) which will be used for the computation of (contracted) GPTs. Note that the background
% conductivity is one. By definition, the n-th diagonal block of A is
% (lambda_n I - K_{D_n}^*), with lambda_n the contrast constant 
% and off diagonal the (m,n)-th block is
% -dn_{m} S_{D_n}.
% What is returned by this function is the blocks of A but the diagonal blocks are replaced by
% K_{D_n}^*. This allows to reuse the result in case of multi-frequency simulation/GPT computation.
%

if ~iscell(D)
    D={D};
end
nbIncls = length(D);

% Check compatibility of inclusions
for m=1:nbIncls
    for n=(m+1):nbIncls
        if ~D{m}.isdisjoint(D{n})
            error('Inclusions must be mutually disjoint.');
        end
    end
end
        
KsdS = cell(nbIncls);

for m=1:nbIncls
    for n=1:nbIncls
        % We use P0-P0 elements, so the stiff matrix is just the kernel matrix
        % Construct Kstar matrix of the n-th inclusion
        if n==m
            KsdS{n,n} = ops.Kstar.make_kernel_matrix(D{n}.points, D{n}.tvec, D{n}.normal, D{n}.avec, ...
                                                        D{n}.sigma);
        else
            KsdS{m,n} = -1*ops.dSLdn.make_kernel_matrix(D{n}.points, D{n}.sigma, D{m}.points, ...
                                                        D{m}.normal);
        end
    end        
end