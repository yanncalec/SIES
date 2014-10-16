%% Class for Single Layer potential operator (Helmholtz equation)
%
classdef SingleLayer_H < ops.Operators
	% Class for single layer potential operator for Helmholtz equation
	
	methods
		function obj = SingleLayer_H(k, D1, type1, step1, D2, type2, step2)
			
			if nargin < 7 % If only one boundary is given
				D2 = D1;
				type2 = type1;
				step2 = step1;
			end
			
			obj = obj@ops.Operators(D1, type1, step1, D2, type2, step2);
			
			if isequal(D1,D2)
				obj.Kmat = ops.SingleLayer_H.make_kernel_matrix(k, D1.points, D1.sigma);
			else
				obj.Kmat = ops.SingleLayer_H.make_kernel_matrix(k, D1.points, D1.sigma, D2.points);
			end
		end
		
	end
	
	methods(Static)
		function K = make_kernel_matrix(k, D, sigma_D, E)
			% Kernel matrix of the single layer potential L^2(\p D) -> L^2(\p E).
			green = @(k,x,y) -1i/4*besselh(0,1,k*sqrt(x.^2+y.^2));
			if nargin == 4
				% Case of two disjoint boundaries, canonical
				% discretization.
				X1 = tools.tensorplus(E(1,:), -D(1,:));
				X2 = tools.tensorplus(E(2,:), -D(2,:));
				K = green(k,X1, X2) .* diag(sigma_D);
			else
				% Case of two identical boundaries. Although the
				% operator is stil well defined (no jumps), the diagonal
				% terms of the kernel matrix are computed in a special way.
				
				N = size(D,2);
				K = zeros(N) ;
				
				for i=1:N
					K(i,1:(i-1)) = sigma_D(1:(i-1)).*green(k,D(1,i)-D(1,1:(i-1)),D(2,i)-D(2,1:(i-1)));
					K(i,(i+1):N) = sigma_D((i+1):N).*green(k,D(1,i)-D(1,(i+1):N),D(2,i)-D(2,(i+1):N));
					% Diagonal terms
					K(i,i) = -1i/4*sigma_D(i)*(2*1i*(-psi(1)+log(sigma_D(i)*k/4) -1)/pi+1);
					
					% The constant \gamma = -psi(1) = 0.5772156649 is the euler's constant%
					%                     Nhan = -1i/4*tools.Helmholtz.diagonal_funct(100000, k, sigma_D(i))
					%                     K(i,i) - Nhan;
					%                     K(i,i) = 2*1i*sigma_D(i)*(log(sigma_D(i)/2)-1)/pi ;
				end
			end
		end
		
		function [ val ] = eval(k, D, F, X)
			% Evaluate the single layer potential S_D[F](X)
			% D: shape
			% F: function defined on D
			% X: points of evaluation, dimension 2X?
			
			G = tools.Helmholtz.Green2D(k, D.points, X);
			val = (reshape(F, 1, []) .* D.sigma) * G;
		end
		
		function [ val ] = eval_grad(k, D, F, X)
			% Gradient of the single layer potential S_D[F](x)
			% D: shape
			% F: the function defined on D
			% X: the points at which to evaluate the gradient
			
			% gradient of the simple layer potential
			[Gx, Gy] = tools.Helmholtz.Green2D_Grad(k, X, D.points);
			
			v1 = Gx * (F(:) .* D.sigma(:));
			v2 = Gy * (F(:) .* D.sigma(:));
			val = transpose([v1 v2]);
		end
	end
	
end

