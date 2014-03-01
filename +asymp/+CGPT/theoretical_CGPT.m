function M = theoretical_CGPT(D, lambda, ord)
% Calculate the CGPT of multiple inclusions D up to a given order.
% The computation is based on the CGPT's integral definition. For a single inclusion:
%
% M_{a,b} = \int_{\partial D} (\lambda - K_D^*)^{-1}[\nu_x \nabla Q(x)](y) P(y) d\delta(y)
%
% with Q(x), P(x) two homogeneous harmonic polynomial corresponding to
% the real or imaginary part of z^n. The CGPT for multiple inclusions is defined similarly, see
% 4.10 of the book Ammari, Kang, 2007.
%
% Inputs:
% D: a list of C2boundary objects
% lambda: constants of contrast
% ord: maximum order of CGPT to be computed
%
% Output:
% M: a matrix of dimension (2*ord) X (2*ord) with the o-o(odd-odd), o-e, e-o,
% e-e(even-even) entries corresponding to respectively CC, CS, SC, and SS matrices.
% 
% Convention:
% M^cc(m,n) = \int_{\p D} Re(z^n) (\lambda I-K_D^*)^-1 [d (Re z^m) dn] ds(x)
%

if ~iscell(D)
    D = {D};
end

nbPoints = D{1}.nbPoints; % all C2boundary objects must have the same value of nbPoints

nbIncls = length(D);
Amat = asymp.CGPT.make_system_matrix(D, lambda);

CC = zeros(ord);
CS = zeros(ord);
SC = zeros(ord);
SS = zeros(ord);

for m=1:ord
    % Right hand vector
    B = zeros(nbPoints, nbIncls);
    for i=1:nbIncls
        dm = m*(D{i}.cpoints).^(m-1); % grad(z^m) 
        toto = D{i}.normal(1,:) .* dm + D{i}.normal(2,:) .* dm * 1i;        
        B(:,i) = toto(:);
    end

    toto = Amat\real(B(:));
    realphim = reshape(toto, nbPoints, nbIncls);
    toto = Amat\imag(B(:));
    imagphim = reshape(toto, nbPoints, nbIncls);
    
    for n=1:ord
        for i=1:nbIncls
            zn = (D{i}.cpoints.^n) .* D{i}.sigma;
            
            CC(m,n) = CC(m,n) + real(zn)*realphim(:,i);
            CS(m,n) = CS(m,n) + imag(zn)*realphim(:,i);
            SC(m,n) = SC(m,n) + real(zn)*imagphim(:,i);
            SS(m,n) = SS(m,n) + imag(zn)*imagphim(:,i);
        end
    end
end

M = asymp.CGPT.cell2mat({CC,CS,SC,SS});

% Symmetrization to improve numerical precision. Note that the complex CGPT matrix is only
% symmetric, not hermitian.
% M = (M+M.')/2;
end

