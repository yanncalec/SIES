%% Class for Double Layer potential
%
% <latex>
% Define:
% \begin{eqnarray}
% D_D[f](x) &=& \int_{\partial D} \frac{\partial}{\partial\nu_y} G(x-y)f(y)ds(y)\\
% &\simeq& \sum_i \langle -\nabla G(x(t_j) - y(t_i)), \nu_y(t_i) \rangle f(y(t_i)) \sigma_D(t_i)
% \end{eqnarray}
% The kernel matrix $K_{ji}$ then reads:
% $$K_{ji} = -\left( \partial_1 G(x(t_j) - y(t_i))\nu^1_y(t_i) + \partial_2
% G(x(t_j) - y(t_i))\nu^2_y(t_i) \right) \sigma_D(t_i)$$
% </latex>
%
classdef DoubleLayer < ops.Operators
	% Double layer potential
	% This operator is not defined for two identical boundaries because of
	% the jump.
	
	methods
		function obj = DoubleLayer(D1, type1, step1, D2, type2, step2)
			if isequal(D1,D2)
				error('Type Error: this operator is not defined for two identical boundaries.');
			end
			
			obj = obj@ops.Operators(D1, type1, step1, D2, type2, step2);
			obj.Kmat = ops.DoubleLayer.make_kernel_matrix(D1.points, D1.normal, D1.sigma, D2.points);
		end
		
		function val = fwd(obj, f)
			if isequal(D1,D2)
				error('Type Error: this operator is not defined for two identical boundaries.');
			else
				val = fwd@ops.Operators(obj, f);
			end
		end
	end
	
	methods(Static)
		function K = make_kernel_matrix(D, normal_D, sigma_D, E)
			% Kernel matrix of the double layer potential L^2(\p D) ->
			% L^2(\p E).
			% Inputs:
			% D, normal_D, sigma_D: boundary, its normal vector and the
			% integration elements.
			% E: boundary.
			
			[Gx, Gy] = tools.Laplacian.Green2D_Grad(E, D);
			K = -(Gx * diag(normal_D(1,:)) + Gy * diag(normal_D(2,:))) * diag(sigma_D);
		end
		
		function [ val ] = eval(D, F, X)
			% Evaluate the double layer potential D_D[F](X)
			% D: shape
			% F: function defined on D
			% X: points of evaluation, dimension 2X?
			
			%             [Gx Gy] = tools.Laplacian.Green2D_Grad(D.points, X);
			%             dGn = [diag(D.normal(1,:)), diag(D.normal(2,:))] * [Gx; Gy];
			% val = (reshape(F, 1, []) .* D.sigma) * dGn;
			
			dGn = tools.Laplacian.Green2D_Dn(X, D.points, D.normal);
			val = dGn * (F(:) .* D.sigma(:)) ;
			val = reshape(val, 1, []);
		end
		
		function [ val ] = eval_grad(D, F, X)
			% Gradient of the double layer potential D_D[F](x)
			% D: shape
			% F: the function defined on D
			% X: the points at which to evaluate the gradient
			
			[H, ~] = tools.Laplacian.Green2D_Hessian(D.points, X);
			vv = -(D.sigma .* reshape(F,1,[])) * (tools.bdiag(D.normal', 2) * H);
			
			% Or equivalently
			% [~, H] = tools.Laplacian.Green2D_hessian(D.points, X);
			% ddGn = [diag(D.normal(1,:)), diag(D.normal(2,:))] * H;
			% M = size(ddGn,2)/2;
			%
			% % perturb the column
			% DD = zeros(size(ddGn));
			% DD(:, 1:2:end) = ddGn(:, 1:M);
			% DD(:, 2:2:end) = ddGn(:, M+1:end);
			%
			% vv = -(D.sigma .* F) * DD;
			
			val = [vv(1:2:end); vv(2:2:end)];
		end
	end
end

