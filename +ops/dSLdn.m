%% Class for $\frac{\partial}{\partial\nu}S_D$
% 
% <latex>
% Define:
% \begin{eqnarray}
% \frac{\partial}{\partial\nu}S_D[f](x) &=& \int_{\partial D} \frac{\partial}{\partial\nu_x} G(x-y)f(y)ds(y)\\
% &\simeq& \sum_i \langle \nabla G(x(t_j) - y(t_i)), \nu_x(t_j) \rangle f(y(t_i)) \sigma_D(t_i)
% \end{eqnarray}
% The kernel matrix $K_{ji}$ then reads:
% $$K_{ji} = \left( \partial_1 G(x(t_j) - y(t_i))\nu^1_x(t_j) + \partial_2
% G(x(t_j) - y(t_i))\nu^2_x(t_j) \right) \sigma_D(t_i)$$
% </latex>
% 

classdef dSLdn < ops.Operators
% Normal derivative of the single layer potential.
%
% This operator is not defined for two identical boundaries because of
% the jump.
    
    methods
        function obj = dSLdn(D1, type1, step1, D2, type2, step2)
            if isequal(D1,D2)
                error('Type Error: this operator is not defined for two identical boundaries.');
            end
            
            obj = obj@ops.Operators(D1, type1, step1, D2, type2, step2);
            obj.Kmat = ops.dSLdn.make_kernel_matrix(D1.points, D1.sigma, D2.points, D2.normal);
        end        
    end
    
    methods(Static)
        function K = make_kernel_matrix(D, sigma_D, E, normal_E)
        % Kernel matrix of the normal derivative of the single layer
        % potential L^2(\p D) -> L^2(\p E).
        % Inputs:
        % D, sigma_D: boundary, and the integration elements.
        % E, normal_E: boundary and its normal vectors.
            
            [Gx, Gy] = tools.Laplacian.Green2D_Grad(E, D);
            K = (diag(normal_E(1,:))*Gx + diag(normal_E(2,:))*Gy) * diag(sigma_D);
        end
        
        function val = eval()
            error('Method not defined because of the jump!');
        end
        
    end
    
end

