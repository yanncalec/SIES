function CGPT = cell2mat(CGPTc)
% CGPT = cell2mat(CGPTc)
% Put the CGPT into a matrix form.  
% Input: 
% CGPTc: CGPT in cell form, with the four elements
% corresponding to {cc, cs, sc, ss} CGPTs.
% Output:
% CGPT: a matrix of dimension (2*ord) X (2*ord) with the
% o-o(odd-odd), o-e, e-o, e-e(even-even) entries corresponding to
% respectively CC, CS, SC, and SS matrices.

CGPT = zeros(2*size(CGPTc{1}));

CGPT(1:2:end, 1:2:end) = CGPTc{1}; % CC
CGPT(1:2:end, 2:2:end) = CGPTc{2}; % CS
CGPT(2:2:end, 1:2:end) = CGPTc{3}; % SC
CGPT(2:2:end, 2:2:end) = CGPTc{4}; % SS
