%% Class for Adjoint of the Neumann-Poincaré operator $\mathcal{K}_D^*$
%
% <latex>
% Define:
% $$ \mathcal{K}_D^*[f](x) = \int_{\partial D} \frac{\langle y-x, \nu_y
% \rangle}{\abs{x-y}^2} f(y) ds(y)$$
% </latex>
%

classdef Kstar < ops.Operators
	% Adjoint of the Neumann-Poincaré operator
	
	methods
		function obj = Kstar(D, type, step)
			obj = obj@ops.Operators(D, type, step, D, type, step);
			
			obj.Kmat = ops.Kstar.make_kernel_matrix(D.points, D.tvec, D.normal, D.avec, D.sigma);
		end
		
	end
	
	methods(Static)
		function Ks = make_kernel_matrix(D, tvec, normal, avec, sigma)
			M = size(D,2);
			Ks = zeros(M, M);
			tvec_norm_square = tvec(1,:).^2 + tvec(2,:).^2;
			
			for j = 1:M
				xdoty = (D(1,j)-D(1,:))*normal(1,j)+(D(2,j)-D(2,:))*normal(2,j);
				norm_xy_square = (D(1,j)-D(1,:)).^2+(D(2,j)-D(2,:)).^2;
				
				Ks(j, 1:j-1) = 1/(2*pi)*xdoty(1:j-1).*(sigma(1:j-1)./norm_xy_square(1:j-1));
				Ks(j, j+1:M) = 1/(2*pi)*xdoty(j+1:M).*(sigma(j+1:M)./norm_xy_square(j+1:M));
			end
			
			for j = 1:M
				%    Ks(j,j) = 1/(2*pi)*(-1)/2*avec(:,j)'*normal(:,j)/tvec_norm_square(j)*sigma(j);
				Ks(j,j) = 1/(2*pi)*(-1)/2*avec(:,j)'*normal(:,j)/tvec_norm_square(j)*sigma(j);
			end
		end
		
		function [ val ] = eval(D, F)
			% Evaluate K_D^*[F]
			% INPUTS:
			% D: shape, a C2boundary object
			% F: function defined on D
			
			Ks = ops.Kstar.make_kernel_matrix(D.points, D.tvec, D.normal, D.avec, ...
				D.sigma);
			val = Ks*F(:);
		end
		
	end
end

