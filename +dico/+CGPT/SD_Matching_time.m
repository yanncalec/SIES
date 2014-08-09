function [err, idx] = SD_Matching_time(I, D, sidx) 
% Dictionary matching algorithm by shape descriptor
% Inputs:
% I: a cell of matrix of shape descriptors to be identified. Each I{m} has dimension (Ntime X scl) and corresponds
% to an independent observation.
% D: a cell of matrix of shape descriptors of dictionary. Each D{n} has
% dimension (Ntime X scl) and corresponds to the n-th shape of the dictionary
% sidx: active scales used for matching, optional
% Outputs:
% err: err(m,:) contains the similarity of the m-th data and the dictionary
% idx: sorted dictionary elements in decreasing order of similarity

if ~iscell(I)
    I = {I}; 
    %error('Input shape descriptors I must be a cell!');
end

M = length(I);
N = length(D);

Ntime = size(I{1}, 1);

if nargin<3
    sidx = 1:size(I{1}, 2);
end
scl = length(sidx);

err=zeros(M, N);
idx=zeros(M, N);

for m=1:M
    for n=1:N
        E = I{m}(:,sidx) - D{n}(:,sidx);
    
        err(m, n) = sqrt(sum(E(:).*E(:)) /scl/Ntime);
    end
    [~, toto] = sort(err(m,:));    
    idx(m,:) = toto;
end

%% Old version: SD is a 3D matrix
% function [err, idx] = SD_Matching_time(I, D, sidx) 
% % Dictionary matching algorithm by shape descriptor
% % Inputs:
% % I: a cell of matrix of shape descriptors to be identified. Each I{m} has dimension (Ntime X scl X 2) and corresponds
% % to an independent observation.
% % D: a cell of matrix of shape descriptors of dictionary. Each D{n} has
% % dimension (Ntime X scl X 2) and corresponds to the n-th shape of the dictionary
% % sidx: active scales used for matching, optional
% % Outputs:
% % err: err(m,:) contains the similarity of the m-th data and the dictionary
% % idx: sorted dictionary elements in decreasing order of similarity
% 
% if ~iscell(I)
%     I = {I}; 
%     %error('Input shape descriptors I must be a cell!');
% end
% 
% M = length(I);
% N = length(D);
% 
% Ntime = size(I{1}, 1);
% 
% if nargin<3
%     sidx = 1:size(I{1}, 2);
% end
% scl = length(sidx);
% 
% err=zeros(M, N);
% idx=zeros(M, N);
% 
% for m=1:M
%     for n=1:N
%         E = I{m}(:,sidx,:) - D{n}(:,sidx,:);
%     
%         err(m, n) = sqrt(sum(E(:).*E(:)) /scl/Ntime);
%     end
%     [~, toto] = sort(err(m,:));    
%     idx(m,:) = toto;
% end
