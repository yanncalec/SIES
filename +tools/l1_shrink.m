function Y = l1_shrink(X, mu, weight)
% l1 reweighted shrinkage.
% Solve min_u mu * \sum_i w_i|u_i| + 1/2 * ||u - v||^2, 
% where ||.|| is the l2 norm and |.| is the abs value. 
% u_i, v_i, w_i are simply scalars, for i=1..N, and 
% |u-v|^2 = \sum_i (u_i - v_i)^2, w_i >= 0
% The solution is given by shrinkage :
% u_i = max(|v_i| - w_i*mu, 0) * v_i / |v_i|
%
% INPUTS:
% X: a vector/matrix 
% mu: 
% weight:
    
    if nargin < 3 || isempty(weight)
        weight = 1;
    end

    Y = zeros(size(X)); 
    V = abs(X) - mu * weight;
    idx = find(V>0);
    Y(idx) = V(idx) .* sign(X(idx));

end