function Y = SXR_op(X,S,R,transp_flag)
% This function implements the linear operator:
% L(M) := S M R^H, (^H is the hermitian of a matrix)
% S and R are matrices of source and receiver
% and the adjoint of L is:
% L^*(D) := S^H D R 

[Ns,Ms]=size(S);
[Nr,Mr]=size(R);

if strcmp(transp_flag,'notransp')
    Y = S * reshape(X, Ms, Mr) * R';
    Y=Y(:);
elseif strcmp(transp_flag,'transp')
    Y = S' * reshape(X, Ns, Nr) * R;
    Y=Y(:);
end
