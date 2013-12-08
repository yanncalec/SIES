function Y = PhaseLift_op(A, X, tflag)
% Operator involved in the PhaseLift technique
% Following the Sparco convention

    [M, N]=size(A);
    % M=round(sqrt(M0));

    % if M*M ~= M0
    %     error('Dimension error: matrix A only apply on square matrix!');
    % end
    
    if tflag==0
        Y={[N, N],[M,1]};
    elseif tflag==1
        X = reshape(X,N,N);
        Y = diag(A*X*A');
    else
        Y = A'*diag(X(:))*A;
        % Y=Y(:);
    end

    % if strcmp(tflag, 'notransp')
end
