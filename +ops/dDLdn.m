%% Class for $\frac{\partial}{\partial\nu}D_D$
%
% <latex>
% Suppose $x\notin \partial D$ and define:
% \begin{equation}
% \frac{\partial}{\partial\nu}D_D[f](x) = \int_{\partial D} \frac{\partial}{\partial\nu_x} \left( \frac{\partial}{\partial\nu_y} G(x-y)\right) f(y)ds(y)
% \end{equation}
% Denote by $D^2G(x)$ the hessian matrix of $G$ evaluated at $x$, then
% $$\frac{\partial}{\partial\nu_x} \left( \frac{\partial}{\partial\nu_y}
% G(x-y)\right) = -\langle \nu_x, D^2G(x-y) \nu_y \rangle,$$
% therefore the integral is approximated by:
% $$-\sum_i \langle \nu_x(t_j), D^2G(x(t_j)-y(t_i)) \nu_y(t_i) \rangle
% f(y(t_i)) \sigma_D(t_i)$$
% The kernel matrix $K_{ji}$ then reads:
% $$K_{ji} = -\langle \nu_x(t_j), D^2G(x(t_j)-y(t_i)) \nu_y(t_i) \rangle \sigma_D(t_i)$$
% If $x\in \partial D$, this operator is \emph{hypersingular} and we have:
% $$\int_{\partial D} \frac{\partial D_D[v_1]}{\partial\nu} v_2 =
% \int_{\partial D} S_D[v_1'] v_2'$$
% </latex>
%

classdef dDLdn < ops.Operators
% Normal derivative of the double layer potential
% 
% This operator is defined through the single layer potential in the hyper singular
% case (D1==D2), and the stiffness matrix exists only for (at least) P1 elements.
    
    properties
        hypersgl=0 % hyper singluar operator: true if D1==D2
    end
    
    methods
        function obj = dDLdn(D1, type1, step1, D2, type2, step2)
            if nargin < 4
                D2 = D1;
                type2 = type1;
                step2 = step1;
            end
            
            obj = obj@ops.Operators(D1, type1, step1, D2, type2, step2);
            
            if ~isequal(D1,D2) % D1!=D2, the operator is well defined
                               % dn_x dn_y G(x,y) = -<n_x, D^2 G(x-y).n_y>
                obj.Kmat = ops.dDLdn.make_kernel_matrix(D1.points, D1.normal, D1.sigma, D2.points, D2.normal);
            else
                % D1==D2, the operator is hyper singular and P0 elements
                % are inappropriate, and we must redefine the boundary
                % element basis Phi and Psi initiated by the constructor of
                % the super class.
                obj.hypersgl = 1;
                
                switch obj.typeBEM1
                  case 'P0'
                    error('Type Error: P0 boundary element is inappropriate for this operator.');
                  case 'P1'
                    %obj.nbBEM1 = floor(D1.nbPoints / step1);
                    obj.Psi = tools.BEM.P1_derivative(D1.nbPoints, step1, 2*pi);
                  otherwise
                    error('Type Error: only P1 elements are supported in the current version.');
                end
                
                switch obj.typeBEM2
                  case 'P0'
                    error('Type Error: P0 boundary element is inappropriate for this operator.');
                  case 'P1'
                    %obj.nbBEM2 = floor(D2.nbPoints / step2);
                    obj.Phi = tools.BEM.P1_derivative(D2.nbPoints, step2, 2*pi);
                  otherwise
                    error('Type Error: only P1 elements are supported in the current version.');
                end
                obj.Psi_t = obj.Psi';
                obj.Phi_t = obj.Phi';
                obj.Kmat = []; % No kernel matrix in this case
            end
        end
        
        function val = fwd(obj, f)
            if isequal(D1,D2)
                error('Type Error: this operator is not defined for two identical boundaries.');
            else
                val = fwd@ops.Operators(obj, f);
            end
        end
        
        function val = apply_stiffmat(obj, f)
            K = obj.get_stiffmat(obj);
            val = K*f;     
        end
        
        function val = get_stiffmat(obj)
            if obj.hypersgl
                Smat = ops.SingleLayer.make_kernel_matrix(obj.D1.points, obj.D1.sigma);
                %val = obj.Phi_t*Smat*diag(obj.D1.sigma)*obj.Psi;
                val = obj.Phi_t * diag(obj.D1.sigma) * (Smat*obj.Psi);
            else
                val = get_stiffmat@ops.Operators(obj);
            end
        end
    end
    
    methods(Static)
        function K = make_kernel_matrix(D, normal_D, sigma_D, E, normal_E)
            [~, hessian] = tools.Laplacian.Green2D_Hessian(E, D);
            H = -[diag(normal_E(1,:)) diag(normal_E(2,:))]*hessian*[diag(normal_D(1,:)); diag(normal_D(2,:))];
            K = H * diag(sigma_D);
        end
        
        function val = eval()
            error('Method not implemented!');
        end
        
    end
end

