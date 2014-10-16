function CGPT_tsr = CGPT_inverse_transform(CGPT, T0, S0, Phi0)
% On the CGPT, apply first the translation -T0, then the scaling 1/S0
% and finally the rotation -Phi0. The output is the transformed CGPT.
%

if iscell(CGPT)
	CGPT=asymp.CGPT.cell2mat(CGPT);
end

[N1,N2]=asymp.CGPT.CGPT2CCGPT(CGPT);
[Z1,Z2]=asymp.CCGPT_inverse_transform(N1,N2,T0,S0,Phi0);

CGPT_tsr=asymp.CGPT.CCGPT2CGPT(Z1,Z2);