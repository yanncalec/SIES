%% Class for Single Layer potential operator
%
% <latex>
% Define:
% $$S_D[f](x) = \int_{\partial D} G(x-y)f(y)ds(y)$$
% which is well defined in the whole space. In case that $x\notin \partial D$, the kernel matrix
% $$K_{ji} = G(x(t_j)-y(t_i)) \sigma_1(t_i),$$
% otherwise the diagonal term of the kernel matrix becomes:
% $$\frac{1}{2\pi} \sigma_1(t_i)\left(\log{\frac{\sigma_1(t_i)}{2}} - 1\right)$$
% </latex>
%

classdef SingleLayer < ops.Operators
% Class for single layer potential operator.
    
    methods
        function obj = SingleLayer(D1, type1, step1, D2, type2, step2)
            if nargin < 4 % If only one boundary is given
                D2 = D1;
                type2 = type1;
                step2 = step1;
            end
            
            obj = obj@ops.Operators(D1, type1, step1, D2, type2, step2);
            
            if isequal(D1,D2)
                obj.Kmat = ops.SingleLayer.make_kernel_matrix(D1.points, D1.sigma);
            else
                obj.Kmat = ops.SingleLayer.make_kernel_matrix(D1.points, D1.sigma, D2.points);
            end
        end
        
    end
    
    methods(Static)
        function K = make_kernel_matrix(D, sigma_D, E)
        % Kernel matrix of the single layer potential L^2(\p D) -> L^2(\p E).
        % Inputs:
        % D, sigma_D: boundary, and the integration elements.
        % E: boundary, optional. If given then E must be disjoint of D.
            
            green = @(x,y) 1/4/pi * log(x.^2+y.^2) ; % Green function
            if nargin == 3
                % Case of two disjoint boundaries, canonical
                % discretization.
                X1 = tools.tensorplus(E(1,:), -D(1,:));
                X2 = tools.tensorplus(E(2,:), -D(2,:));
                K = green(X1, X2) * diag(sigma_D);
            else
                % Case of two identical boundaries. Although the
                % operator is stil well defined (no jumps), the diagonal
                % terms of the kernel matrix are computed in a special way.
                
                N = size(D,2);
                K = zeros(N) ;
                
                for i=1:N
                    K(i,1:i-1) = sigma_D(1:i-1).*green(D(1,i)-D(1,1:i-1),D(2,i)-D(2,1:i-1));
                    K(i,i+1:N) = sigma_D(i+1:N).*green(D(1,i)-D(1,i+1:N),D(2,i)-D(2,i+1:N));
                    % Diagonal terms
                    % K(i,i) = sigma(i)*(log(sigma(i))-1)/2/pi ; % Thomas' original version
                    K(i,i) = sigma_D(i)*(log(sigma_D(i)/2)-1)/2/pi ; % Han's version
                end
            end
        end
        
        function [ val ] = eval(D, F, X)
        % Evaluate the single layer potential S_D[F](X)
        % INPUTS:
        % D: shape, a C2boundary object
        % F: function defined on D
        % X: points of evaluation of dimension 2X?. Must not contain any points of D.
            
            G = tools.Laplacian.Green2D(D.points, X);           
            val = (reshape(F, 1, []) .* D.sigma) * G;
        end
        
        function [ val ] = eval_grad(D, F, X)
        % Gradient of the single layer potential S_D[F](x)
        % INPUTS:
        % D: shape, a C2boundary object
        % F: function defined on D
        % X: points of evaluation of dimension 2X?. Must not contain any points of D.

        % gradient of the simple layer potential
            [Gx, Gy] = tools.Laplacian.Green2D_Grad(X, D.points);
            
            v1 = Gx * (F(:) .* D.sigma(:));
            v2 = Gy * (F(:) .* D.sigma(:));
            val = transpose([v1 v2]);
        end
    end
    
end

