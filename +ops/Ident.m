classdef Ident < ops.Operators
	% Identity operator
	
	methods
		function obj = Ident(D, type1, step1, type2, step2)
			if nargin < 4
				type2 = type1;
				step2 = step1;
			end
			obj = obj@ops.Operators(D, type1, step1, D, type2, step2);
			obj.Kmat = eye(D.nbPoints);
		end
	end
	
	methods(Static)
		function K = make_kernel_matrix(M)
			K = eye(M);
		end
		
		function val = eval(F)
			val = F;
		end
	end
end

