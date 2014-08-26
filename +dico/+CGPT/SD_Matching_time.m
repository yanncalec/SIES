function [err, idx] = SD_Matching_time(I, D, Sidx) 
% Dictionary matching algorithm by shape descriptor
% Inputs:
% I: a cell of matrix of shape descriptors to be identified. Each I{m} has dimension (Ntime X scl) and corresponds
% to an independent observation.
% D: a cell of matrix of shape descriptors of dictionary. Each D{n} has
% dimension (Ntime X scl) and corresponds to the n-th shape of the dictionary
% Sidx: active scales used for matching, optional
% Outputs:
% err: err(m,:,s) contains the similarity of the m-th data and the dictionary. The third
% dimension is the comparaison using the first s scales in Sidx
% idx: sorted dictionary elements in decreasing order of similarity

if ~iscell(I)
    I = {I}; 
    %error('Input shape descriptors I must be a cell!');
end

M = length(I);
N = length(D);

Ntime = size(I{1}, 1);

if nargin<3
    Sidx = 1:size(I{1}, 2);
end

nbScl = length(Sidx);

err=zeros(M, N, nbScl);
idx=zeros(M, N, nbScl);

for m=1:M
    for n=1:N
        E = I{m}(:,Sidx,:) - D{n}(:,Sidx,:);
        Esq = sum(E.*E, 1);

        for s=1:nbScl
            err(m, n, s) = sqrt(sum(Esq(1:s))/s/Ntime);
            [~, idx(m,:,s)] = sort(err(m,:,s));
        end    
    end
end

