function Y = SXR_op_list(X,S,R,transp_flag)
% Y = SXR_op_list(X,S,R,transp_flag)
% Implementation of the following operator
% L(X)_{sr} = \sum_{mn} S_{sm} X_{mn} (R(s)^H)_{nr}
% with R a list of matrices (R(s) is the s-th matrix),
% and its adjoint:
% L^*(Z)_{mn} = \sum_{sr} S^H_{ms} Z_{sr} R(s)_{rn}
%
% This operator corresponds to the acquisition using source dependant
% (non-fixed) receivers, eg, the ElectricFish problem.

[Ns,Ms]=size(S);
[Nr,Mr]=size(R{1});

if strcmp(transp_flag,'notransp')
    Xm = reshape(X, Ms, Mr);
    Y = zeros(Ns, Nr);
    for s=1:Ns
        Y(s, :) = S(s,:) * Xm * R{s}';
    end
    Y=Y(:);
elseif strcmp(transp_flag,'transp')
    Z = reshape(X, Ns, Nr); 
    W = zeros(Ns, Mr);
    for s=1:Ns
        W(s, :) = Z(s,:) * R{s};
    end
    Y = S' * W;
    Y=Y(:);
end
