function [N1, N2] = CGPT2CCGPT(CGPT)
% Make the Complex CGPT from the CGPT. 
% Input:
% CGPT: the standard CGPT in matrix or cell form
% Outputs:
% N1, N2: the Complex CGPT such that N1=N1^T, N2=N2^H

if ~iscell(CGPT)
    CGPT = asymp.CGPT.mat2cell(CGPT);
end

N1 = CGPT{1} - CGPT{4} + 1i * (CGPT{2} + CGPT{3});
N2 = CGPT{1} + CGPT{4} + 1i * (CGPT{2} - CGPT{3});

% % Symmetrization
% N1 = (N1 + N1.')/2; % N1 is symmetric
% N2 = (N2 + N2')/2; % N2 is hermitian
