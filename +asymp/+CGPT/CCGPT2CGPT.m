function CGPT = CCGPT2CGPT(N1,N2)
% Convert the Complex CGPT matrix N1 and N2 to CGPT
% Input:
% CCGPT: the Complex CGPT
% Outputs:
% CGPT: the standard CGPT in matrix form
%
% WARNING: 
% the conversion to CGPT is valide only for real conductivity case, ie, lambda is a real number. 


CGPT={};

CGPT{1}=real(N1+N2)/2; %CC
CGPT{2}=imag(N1+N2)/2; %CS
CGPT{3}=imag(N1-N2)/2; %SC
CGPT{4}=real(N2-N1)/2; %SS

CGPT=asymp.CGPT.cell2mat(CGPT);

