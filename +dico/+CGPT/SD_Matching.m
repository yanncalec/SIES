function [err, idx] = SD_Matching(I1, I2, D1, D2, cord, method) 
% Dictionary matching algorithm by shape descriptor
% Inputs:
% I1, I2: shape descriptor of data
% D1, D2: shape descriptors of dictionary in cell format
% cord: maximum order for comparison
% method: method for measuring the similarity (see the code)
% Outputs:
% err: cell of similarity. err{k} contains the similarity
% of the first k orders (k=1..cord) descriptors between data and the dictionary
% idx: cell of sorted dictionary elements in decreasing order of similarity

if nargin<6
    method=0;
end

for ord=1:cord
    err{ord}=zeros(length(D1),1);
end

for n=1:length(D1)
    E1 = I1(1:cord, 1:cord)-D1{n}(1:cord, 1:cord); 
    E2 = I2(1:cord, 1:cord)-D2{n}(1:cord, 1:cord);

    for ord=1:cord
        F1 = E1(1:ord, 1:ord);
        F2 = E2(1:ord, 1:ord);

        % Different methods for measuring the similarity
        if method==0
            err{ord}(n) = sqrt(norm(F1, 'fro')^2+norm(F2, 'fro')^2);
        else
            err{ord}(n) = max(norm(F1, 'fro'), norm(F2, 'fro'));
            % err{ord}(n) = min(norm(F1, 'fro'), norm(F2, 'fro'));
        end
        %err(n) = sqrt(norm(F1, 'fro')^2+norm(F2, 'fro')^2) / sqrt(norm(D1, 'fro')^2+norm(D2, 'fro')^2);
        %err(n) = norm(F2, 'fro');
    end
end

% For each order from 1 to cord, the dictionary element is sorted
for ord=1:cord
    % The first one in idx{ord} is the identified element index
    [~,idx{ord}]=sort(err{ord});
end
