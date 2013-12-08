function Y = SXR_op_sym(X,S,R,transp_flag)
% This function implements the linear operator:
% L(M) := S (M+M^T) R^H
% and the adjoint of L is:
% L^H(D) := S^H D R + R^T D^T conj(S)

% The dimension of data matrix
[Ms,Ns]=size(S);
[Mr,Nr]=size(R);

if strcmp(transp_flag,'notransp')
    Xm = reshape(X, Ns, Nr); 
    Y = S * (Xm + Xm.') * R.';
    Y=Y(:);
elseif strcmp(transp_flag,'transp')
    Xm = reshape(X, Ms, Mr);
    Y = S' * Xm * R + R.' * Xm.' * conj(S);
    Y=Y(:);
end
