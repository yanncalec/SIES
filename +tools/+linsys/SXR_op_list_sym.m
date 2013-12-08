function Y = SXR_op_list_sym(X,S,R,transp_flag)
% symmetric version of SXR_op_list

[Ns,Ms]=size(S);
[Nr,Mr]=size(R{1});

if strcmp(transp_flag,'notransp')
    Xm = reshape(X, Ms, Mr); 
    Xm = Xm + Xm.';
    Y = zeros(Ns, Nr);
    for s=1:Ns
        Y(s, :) = S(s,:) * Xm * R{s}.';
    end
    Y=Y(:);
elseif strcmp(transp_flag,'transp')
    Z = reshape(X, Ns, Nr); 
    Zt = Z.';
    W1 = zeros(Ns, Mr);
    W2 = zeros(Mr, Ns);
    for s=1:Ns
        W1(s, :) = Z(s,:) * R{s};
        W2(:, s) = R{s}.' * Zt(:,s);
    end
    Y = S' * W1 + W2 * conj(S);
    Y=Y(:);
end
