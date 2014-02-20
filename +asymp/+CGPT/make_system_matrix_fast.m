function Amat = make_system_matrix_fast(KsdS, lambda)
% Construct the matrix A in the system A[phi]=b by reusing the block matrices
% constructed by the function make_block_matrix.

nbIncls = size(KsdS,1);

if iscell(lambda)
    lambda = cell2mat(lambda);
end

if length(lambda) < nbIncls
    error('Value of lambda must be specified for each inclusion.');
end

Amat = KsdS;

for n = 1:nbIncls
    Amat{n,n} = lambda(n)*eye(size(KsdS{n,n})) + KsdS{n,n}; % It is a plus here by construction of KsdS (see make_block_matrix function)
end

Amat = cell2mat(Amat);
end