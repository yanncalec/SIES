function Amat = make_system_matrix(D, lambda)
% Construct some blocks of the matrix A in the system A[phi]=b

if ~iscell(D)
    D={D};
end
nbIncls = length(D);

if length(lambda) < nbIncls
    error('Value of lambda must be specified for each inclusion.');
end

% Check compatibility of inclusions
for m=1:nbIncls
    for n=(m+1):nbIncls
        if ~D{m}.isdisjoint(D{n})
            error('Inclusions must be mutually disjoint.');
        end
    end
end
        
A = cell(nbIncls);

for m=1:nbIncls
    for n=1:nbIncls
        % We use P0-P0 elements, so the stiff matrix is just the kernel matrix
        % Construct Kstar matrix of the n-th inclusion
        if m==n
            toto = ops.Kstar.make_kernel_matrix(D{n}.points, D{n}.tvec, D{n}.normal, D{n}.avec, ...
                                                D{n}.sigma);
            A{n,n} = lambda(n)*eye(size(toto)) - toto;
        else
            A{m,n} = -1*ops.dSLdn.make_kernel_matrix(D{n}.points, D{n}.sigma, D{m}.points, ...
                                                        D{m}.normal);
        end
    end        
end

Amat = cell2mat(A);
