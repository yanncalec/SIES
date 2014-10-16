%% Class for boundary operators
%
% For documentation see the latex file Operators.tex

%%
classdef Operators < handle
	% Abstract class for linear (integral) operators defined on two
	% boundaries. This class provides interfaces for the stiffness matrix
	% and the application of the operator.
	
	properties(SetAccess = protected)
		D1 % boundary of the domain of the operator
		typeBEM1 % type of boundary elements for D1
		nbBEM1 % number of boundary elements for D1 (GET)
		stepBEM1 % step length for D1 boundary elements
		Psi % Basis for D1
		Psi_t % adjoint of Psi
		
		D2 % boundary of the image of the operator
		typeBEM2 % type of boundary elements for D2
		nbBEM2 % number of boundary elements for D2 (GET)
		stepBEM2 % step length for D2 boundary elements
		Phi % Basis for D2
		Phi_t % adjoint of Phi
		
		Kmat % kernel matrix
		stiffmat% stiffness matrix representation of the operator. stiffmat is the same as Kmat under the P0-P0 basis.
	end
	
	methods(Abstract, Static)
		K = make_kernel_matrix() % Compute the kernel matrix of the operator
		val = eval() % evaluate the operator
	end
	
	methods
		function obj = Operators(D1, type1, step1, D2, type2, step2)
			obj.D1 = D1;
			obj.typeBEM1 = type1;
			obj.D2 = D2;
			obj.typeBEM2 = type2;
			
			switch obj.typeBEM1
				case 'P0'
					% If P0 basis then the dimension of the approximation space equals to the number of
					% boundary points
					
					% obj.Psi = [];
					obj.Psi = eye(D1.nbPoints);
					obj.stepBEM1 = 1;
				case 'P1'
					if mod(D1.nbPoints, step1)
						error('Invalid sampling step for D1');
					end
					obj.stepBEM1 = step1;
					obj.Psi = tools.BEM.P1_basis(D1.nbPoints, step1);
				otherwise
					error('Type Error: only P0 and P1 elements are supported in the current version.');
			end
			
			switch obj.typeBEM2
				case 'P0'
					% obj.Phi = [];
					obj.Phi = eye(D2.nbPoints);
					obj.stepBEM2 = 1;
				case 'P1'
					if mod(D2.nbPoints, step2)
						error('Invalid sampling step for D2');
					end
					obj.stepBEM2 = step2;
					obj.Phi = tools.BEM.P1_basis(D2.nbPoints, step2);
				otherwise
					error('Type Error: only P0 and P1 elements are supported in the current version.');
			end
			
			obj.Psi_t = obj.Psi';
			obj.Phi_t = obj.Phi';
		end
		
		function val = get.nbBEM1(obj)
			val = floor(obj.D1.nbPoints / obj.stepBEM1);
		end
		
		function val = get.nbBEM2(obj)
			val = floor(obj.D2.nbPoints / obj.stepBEM2);
		end
		
		function val = fwd(obj, f)
			% Apply the operator $A$ on a $L^2(\partial D_1)$ function $f$
			% represented by its synthesis coefficients wrt the BEM basis
			% $\Psi$. The result is the samples of $Af$ on $\partial D_2$.
			% This is different with the application of stiffness matrix
			% which yields the analysis coefficient of $Af$ wrt the
			% basis on $\partial D_2$.
			val = obj.Kmat * obj.Psi * f;
			
			% if strcmp(obj.typeBEM1, 'P0')
			%     val = obj.Kmat * f ;
			% elseif  ~strcmp(obj.typeBEM1, 'P0')
			%     val = obj.Kmat * obj.Psi * f;
			% end
		end
		
		function val = get.stiffmat(obj)
			val = obj.get_stiffmat();
		end
		
		function val = get_stiffmat(obj)
			% Export the stiffness matrix
			sigma2 = diag(obj.D2.sigma);
			
			if strcmp(obj.typeBEM1, 'P0') && strcmp(obj.typeBEM2, 'P0')
				val = obj.Kmat;
			elseif strcmp(obj.typeBEM1, 'P0') && ~strcmp(obj.typeBEM2, 'P0')
				val = obj.Phi_t * sigma2 * obj.Kmat;
			elseif ~strcmp(obj.typeBEM1, 'P0') && strcmp(obj.typeBEM2, 'P0')
				val = obj.Kmat * obj.Psi;
			elseif ~strcmp(obj.typeBEM1, 'P0') && ~strcmp(obj.typeBEM2, 'P0')
				val = obj.Phi_t * sigma2 * obj.Kmat * obj.Psi;
			else
				error('Not implemented');
			end
			
		end
		
		% function val = apply_stiffmat(obj, f)
		% % Apply the stiffness matrix on a vector without explicitly
		% % construct the matrix
		%     sigma2 = obj.D2.sigma(:);
		
		%     if strcmp(obj.typeBEM1, 'P0') && strcmp(obj.typeBEM2, 'P0')
		%         val = sigma2 .* (obj.Kmat * f) ;
		%     elseif strcmp(obj.typeBEM1, 'P0') && ~strcmp(obj.typeBEM2, 'P0')
		%         val = obj.Phi_t * (sigma2 .* (obj.Kmat * f));
		%     elseif ~strcmp(obj.typeBEM1, 'P0') && strcmp(obj.typeBEM2, 'P0')
		%         val = sigma2 .* (obj.Kmat * (obj.Psi * f));
		%     elseif ~strcmp(obj.typeBEM1, 'P0') && ~strcmp(obj.typeBEM2, 'P0')
		%         val = obj.Phi_t * (sigma2 .* (obj.Kmat * (obj.Psi * f)));
		%     else
		%         error('Not implemented');
		%     end
		% end
		
	end
end
