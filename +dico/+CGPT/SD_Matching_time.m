function [err, idx] = SD_Matching_time(I, D, method) 
% Dictionary matching algorithm by shape descriptor
% Inputs:
% I: a matrix of shape descriptors (data) to be identified, each row corresponds
% to an independent observation
% D: shape descriptors of dictionary in a matrix format with D(n,:) the
% shape descriptor corresponding to the n-th shape
% Outputs:
% err: err(m,:) contains the similarity of the m-th data and the dictionary
% idx: sorted dictionary elements in decreasing order of similarity

if nargin<3
    method = 0;
end

N = size(D,1);
M = size(I,1);

err=zeros(M, N);
idx=zeros(M, N);

for m=1:M
    for n=1:N
        E = I(m,:) - D(n, :);
    
        % Different methods for measuring the similarity
        if method==0
            err(m, n) = norm(E, 'fro');
        else
            0;
        end
    end
    [~, toto] = sort(err(m,:));    
    idx(m,:) = toto;
end

