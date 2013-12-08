classdef dSLdn_H < ops.Operators
% Normal derivative of the single layer potential for Helmholtz equation.
    
    methods
        function obj = dSLdn_H(k, D1, type1, step1, D2, type2, step2)
            if isequal(D1,D2)
                error('Type Error: this operator is not defined for two identical boundaries.');
            end
            
            obj = obj@ops.Operators(D1, type1, step1, D2, type2, step2);
            % obj.k = k;
            obj.Kmat = ops.dSLdn_H.make_kernel_matrix(k, D1.points, D1.sigma, D2.points, D2.normal);
        end
    end
    
    methods(Static)
        function K = make_kernel_matrix(k, D, sigma_D, E, normal_E)
        % Kernel matrix of the normal derivative of the single layer
        % potential L^2(\p D) -> L^2(\p E).
        % Inputs:
        % D, sigma_D: boundary, and the integration elements.
        % E, normal_E: boundary and its normal vectors.
            
            [Gx, Gy] = tools.Helmholtz.Green2D_Grad(k,E,D);
            K = (diag(normal_E(1,:))*Gx + diag(normal_E(2,:))*Gy) * diag(sigma_D);
        end
        
        function val = eval()
            error('Method not implemented!');
        end
        
    end
    
end

