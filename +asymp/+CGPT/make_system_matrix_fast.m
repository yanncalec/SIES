function [Amat, Acell] = make_system_matrix_fast(KsdS, lambda)
% Construct the matrix A in the system A[phi]=b by reusing the block matrices
% constructed by the function make_block_matrix.
% Inputs:
% -KsdS: a cell of block matrix, returned by the function make_block_matrix
% -lambda: contrast constant of each inclusion, an array

nbIncls = size(KsdS,1);

if iscell(lambda)
	lambda = cell2mat(lambda);
end

if length(lambda) < nbIncls
	error('Value of lambda must be specified for each inclusion.');
end

Acell = KsdS;

for n = 1:nbIncls
	Acell{n,n} = lambda(n)*eye(size(KsdS{n,n})) + KsdS{n,n}; % It is a plus here by construction of KsdS (see make_block_matrix function)
end

Amat = cell2mat(Acell);
end