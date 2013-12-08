function CGPTc = mat2cell(CGPT)
% CGPTc = mat2cell(CGPT)
% Put the CGPT into a cell form.  
% Inputs:
% CGPT: the CGPT matrix defined blockwisely
% Outputs: 
% CGPTc: the CGPT in cell form, with the four elements
% corresponding to {cc, cs, sc, ss} CGPTs.

CGPTc{1} = CGPT(1:2:end, 1:2:end); % CC
CGPTc{2} = CGPT(1:2:end, 2:2:end); % CS
CGPTc{3} = CGPT(2:2:end, 1:2:end); % SC
CGPTc{4} = CGPT(2:2:end, 2:2:end); % SS
