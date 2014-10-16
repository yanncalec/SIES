function [err, idx] = PT_Matching(PT, Dico)
% Dictionary matching algorithm by polarization tensor
% Inputs:
% -PT: 2X2 matrix of polarization tensor, data
% -Dico: a dictionary of polarization tensors
% Outputs:
% -err: similarity (in Euclidean norm) between data and the dictionary
% -idx: sorted dictionary elements in decreasing order of similarity

err=zeros(length(Dico),1);

% PT is  invariant to translation
S1 = svd(PT); % Singular values of PT, invariant to rotation
vv = S1(1) / S1(2); % ratio between singular values, invariant to scaling

for n=1:length(Dico)
	S0 = svd(Dico{n});
	vn = S0(1) / S0(2);
	err(n) = norm(vv-vn);
end
% err=zeros(length(Dico),1);
% vv = sort(eig(PT/trace(PT)));

% for n=1:length(Dico)
%     vn = sort(eig(Dico{n}/trace(Dico{n})));
%     err(n) = norm(vv-vn);
% end

[~,idx]=sort(err);
