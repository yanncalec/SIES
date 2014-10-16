%% Class for the operator $(\lambda I - \mathcal{K}_D^*)^{-1}$
%
% <latex>
% This operator is invertible on $L^2(\partial D)$ for $\lambda>1/2$ or
% $\lambda\leq -1/2$, and $L^2(\partial D)$ for $|\lambda|\geq 1/2$.
% </latex>


classdef LmKstarinv < ops.Operators
	methods
		function obj = LmKstarinv(lambda, D, type, step)
			obj = obj@ops.Operators(D, type, step, D, type, step);
			
			obj.Kmat = ops.Kstar.make_kernel_matrix(lambda, D.points, D.tvec, D.normal, D.avec, ...
				D.sigma);
		end
	end
	
	methods(Static)
		function Kmat = make_kernel_matrix(lambda, D, tvec, normal, avec, sigma)
			if abs(lambda)<1/2
				error('The operator is not defined for this value of lambda!');
			end
			
			Ks = ops.Kstar.make_kernel_matrix(D, tvec, normal, avec, sigma);
			A = lambda * eye(size(Ks)) - Ks;
			Kmat = inv(A);
		end
		
		function val = eval()
			error('Method not implemented!');
		end
		
	end
end

