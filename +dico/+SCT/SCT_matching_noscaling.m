function [err, idx] = SCT_matching_noscaling(S, D)
% Dictionary matching algorithm by shape descriptor of SCT without scaling
% INPUTS:
% S: shape-descriptor extracted from data, a 2D array.
% D: the dictionary of pre-computed shape descriptor (inv. to translation and rotation), a
% cell. D{n} is the shape descriptor of the n-th element, a 2d array.
% OUTPUTS:
% err: euclidean distance between S and each element of D
% idx: sorted dictionary elements in decreasing order of similarity

[M,N] = size(S);
[M1,N1] = size(D{1});

if (M~=M1 || N~=N1) 
    error('Dimension error!');
end

err = zeros(1, length(D));

for n=1:length(D)
    err(n) = norm(S-D{n},'fro');
end

[~,idx]=sort(err);
