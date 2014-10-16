classdef Ellipse < shape.C2boundary
	% Class for ellipse shape
	
	properties(SetAccess = protected)
		axis_a = 1; % length of the semi-major axis
		axis_b = 1; % length of the semi-minor axis
		phi = 0; % orientation of the ellipse
	end
	
	methods
		function obj = Ellipse(a, b, nbPoints)
			% a,b: the length of semi-major and semi-minor axis
			
			if a<b
				error(['Value error: the semi-major axis must be longer than the ' ...
					'semi-minor one.']);
			end
			
			com = [0,0]';
			
			% Boundary of the ellipse
			theta = 2*pi*(0:nbPoints-1)/nbPoints;
			points = repmat(com, 1, nbPoints) + [a*cos(theta) ; b*sin(theta)] ;
			
			% Tangent vector of boundary: first order derivative of D(theta)
			tvec = [-a*sin(theta) ; b*cos(theta)] ;
			
			% Acceleration vector: second order derivative of D(theta)
			avec = [-a*cos(theta); -b*sin(theta)] ;
			
			% Outward normal vector of boundary, orthogonal to tvec
			normal = [[0 1];[-1 0]]*tvec ;
			tt = sqrt(normal(1,:).^2+normal(2,:).^2);
			normal = normal ./ repmat(tt,2,1);
			
			if a == b
				name_str = 'Circle';
			else
				name_str = 'Ellipse';
			end
			
			obj = obj@shape.C2boundary(points, tvec, avec, normal, com, name_str);
			
			obj.axis_a = a;
			obj.axis_b = b;
			obj.phi = 0;
		end
		
		function obj = mtimes(obj, s)
			% Overload of the operator *
			obj = mtimes@shape.C2boundary(obj, s);
			obj.axis_a = obj.axis_a * s;
			obj.axis_b = obj.axis_b * s;
		end
		
		function obj = lt(obj, phi)
			% Redefine the < operator as the rotation of the boundary
			obj = lt@shape.C2boundary(obj,phi);
			obj.phi = obj.phi+phi;
		end
		
		function CGPT = CGPT_theo(obj, ord, freq)
			% Theoretical value of the CGPT by explicite formula, well
			% defined only for ellipse centered at 0.
			if nargin < 3
				freq = 0;
			end
			if norm(obj.center_of_mass) > 0
				error('The explicite formula can be used only for the ellipse centered at the origin');
			end
			
			if obj.axis_a == obj.axis_b % in case of a disk
				M = GPT.disktensor(ord, obj.axis_a/2, obj.kappa(freq)); % axis_a is the diameter
			else
				M = GPT.ellipsetensor(ord, obj.axis_a/2, obj.axis_b/2, obj.kappa(freq));
			end
			% M1 = defs.GPT.CGPT_transform(M, obj.center_of_mass, 1, obj.phi);
			CGPT = GPT.CGPT(M);
		end
	end
	
	% methods(Access = protected)
	%     function val = get_center_of_mass(obj)
	%         val = obj.center_of_mass;
	%     end
	% end
end
