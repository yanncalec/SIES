function Y = add_white_noise(X, nlvl, mode, rowmajor)
% Y = add_white_noise(X, nlvl, mode, rowmajor)
%
% Add white noise to data (independantly column by column)
% INPUTS:
% X: input data matrix
% nlvl: noise level
% mode: if mode~=0 then each column of X is treated independantly
% rowmajor: if true then take transpose of X before adding noise (noise is added row by row in lieu of column by
% column by default)
%    
% OUTPUT:
% Y: noisy data
        
    if nargin < 4
       rowmajor = 0;
    end    
    if nargin < 3
        mode = 0;
    end
        
    if ~iscell(X)
        if rowmajor
            Y = transpose(add_white_noise_mat(transpose(X), nlvl, mode));
        else
            Y = add_white_noise_mat(X, nlvl, mode);
        end        
    else
        for n=1:length(X)
            if rowmajor
                Y{n} = transpose(add_white_noise_mat(transpose(X{n}), nlvl, mode));
            else
                Y{n} = add_white_noise_mat(X{n}, nlvl, mode);
            end
        end
    end                            
end

function Y = add_white_noise_mat(X, nlvl, mode)
    [M, N] = size(X);
    
    if mode
        Y = zeros(M,N);
        for c = 1:N
            % Method 1
            % epsilon = mean(abs(X(:,c))) * nlvl;

            % Method 2
            t0 = norm(X(:,c)) / sqrt(M);
            epsilon = t0 * nlvl;
            
            Y(:, c) = X(:,c) + randn(M,1)*epsilon;            
        end
    else
        % Method 1
        % epsilon = mean(abs(X(:))) * nlvl;
        
        % Method 2
        t0 = norm(X,'fro') / sqrt(M*N);
        epsilon = t0 * nlvl;

        Y = X + randn(M,N)*epsilon;
    end
end

% function Y = add_white_noise_mat_complex(X, nlvl, mode)
% %epsilon = (max(X(:)) - min(X(:))) * nlvl;
%     if mode
%         Y = zeros(size(X));
%         for c = 1:size(X,2)
%             epsilon = mean(abs(X(:,c))) * nlvl;
%             RY = real(X(:,c)) + randn(size(X,1),1)*epsilon;            
%             IY = imag(X(:,c)) + randn(size(X,1),1)*epsilon;
%             Y(:, c) = RY + 1i * IY;
%         end
%     else
%         epsilon = mean(abs(X(:))) * nlvl;
%         RY = real(X) + randn(size(X))*epsilon;
%         IY = imag(X) + randn(size(X))*epsilon;
%         Y = RY + 1i * IY;
%     end
% end
