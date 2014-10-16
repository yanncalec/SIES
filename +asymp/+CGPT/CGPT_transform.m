function CGPT_tsr = CGPT_transform(CGPT, T0, S0, Phi0)
% On the CGPT matrix, apply first the rotation Phi0, then the scaling S0,
% finally the translation T0. The output is the transformed CGPT.
%

if iscell(CGPT)
	CGPT=asymp.CGPT.cell2mat(CGPT);
end

[N1,N2]=asymp.CGPT.CGPT2CCGPT(CGPT);
[Z1,Z2]=asymp.CGPT.CCGPT_transform(N1,N2,T0,S0,Phi0);

CGPT_tsr=asymp.CGPT.CCGPT2CGPT(Z1,Z2);