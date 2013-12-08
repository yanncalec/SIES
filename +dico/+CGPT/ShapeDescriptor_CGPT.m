function [I1, I2, S1, S2, T1, T2] = ShapeDescriptor_CGPT(CGPT)
% Make invariant shape descriptor from CGPT
% Inputs:
% N1, N2: the complex CGPT of shape
% Outputs:
% I1, I2: the invariant shape descriptor based on CGPT
% S1, S2: scaling and translation invariant descriptor
% T1, T2: translation invariant descriptor 

    [N1, N2] = asymp.CGPT.CGPT2CCGPT(CGPT);
    [M,N] = size(N1);

    if M~=N || size(N2,1)~=M || size(N2,2)~=N
        error('Square matrices N1 and N2 must have the same size!');
    end

    [T1, T2] = asymp.CGPT.CCGPT_inverse_transform(N1, N2, N2(1,2)/N2(1,1)/2, 1, 0);

    % 1st scaling invariance: instable at high order
    %D = diag(abs(T2(1,1)).^(-(1:M)/2));

    % 2st scaling invariance: stable at high order but the diagonal
    % elements of N2 are identically one.
    D = diag(1./sqrt(abs(diag(T2))));

    S1 = D * T1 * D;
    S2 = D * T2 * D;

    I1 = abs(S1);
    I2 = abs(S2);
end
