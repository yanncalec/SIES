function [err, idx, PT_inv, Dico_inv] = MF_PT_Matching(PT, Dico)
% Multi-frequency dictionary matching using the Polarization Tensors.
% INPUTS:
% -PT : a list of tensors to be matched, PT{n} is a 2-by-2 PT matrix of the n-th frequency
% -Dico : a list of dictionaries, Dico{n}{m} is a 2-by-2 PT matrix of the m-th shape at the n-th frequency
% OUTPUTS:
% -err: similarity (in Euclidean norm) between data and the dictionary
% -idx: sorted dictionary elements in decreasing order of similarity
% -PT_inv, Dico_inv: invariants calculated from PT and Dico
%
% The frequency indexes (and their value) in PT and Dico should be matched, and the frequency should
% be in an increasing order. This algorithm exploits the multi-frequency information. It creats the
% invariance by normalizing the singular value of each frequency against the one of the largest
% frequency. This allows to distinguish rotational symmetric objects which is impossible with a
% single frequency.

Nfreq = length(PT);
NDico = length(Dico{1});

S = zeros(2,Nfreq);
for n=1:Nfreq
	S(:, n) = svd(PT{n});
end
% % Normalization by the maximum frequency
% PT_inv = inv(diag(S(:, end))) * S;
% % Normalization by the mean value
toto = mean(S,2);
PT_inv = inv(diag(toto)) * S;

for m=1:NDico
	S0{m} = zeros(2,Nfreq);
	for n=1:Nfreq
		S0{m}(:, n) = svd(Dico{n}{m});
	end
	% % Normalization by the maximum frequency
	% Dico_inv{m} = inv(diag(S0{m}(:, end))) * S0{m};
	% % Normalization by the mean value
	toto = mean(S0{m},2);
	Dico_inv{m} = inv(diag(toto)) * S0{m};
end

for m=1:NDico
	toto = PT_inv - Dico_inv{m};
	%err{m} = sqrt(toto(1,:).^2 + toto(2,:).^2);
	err(m) = sqrt(norm(toto, 'fro')^2 / Nfreq);
end

[~ , idx] = sort(err);
end
