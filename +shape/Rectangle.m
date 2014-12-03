classdef Rectangle < shape.C2boundary
	% Class for rectangular shape
	properties(SetAccess = protected)
		width
		height
	end
	
	methods
		function obj = Rectangle(a,b,nbPoints,dspl)
			% This function creates a structure representing an anomaly which has the
			% shape of a rectangle.
			% INPUT : - a = height,b = width
			%         - nPoints = number of discretization points
			% dspl: down-sampling factor for smoothing the corners
			
			if nargin<4
				dspl = 10;
			end
			
			t1 = b/(a+b)/2; t2 = a/(a+b)/2; t3 = t1; t4 = t2;
			n1 = floor(t1*nbPoints); n2 = floor(t2*nbPoints);
			n3 = floor(t3*nbPoints); n4 = nbPoints-n1-n2-n3;
			
			A = [-b ; a]/2;
			B = [-b; -a]/2;
			C = [b ; -a]/2;
			D = [b ; a]/2;
			
			t = (0:n1-1)/n1;
			AB = repmat(A,1,n1) + repmat(B-A,1,n1) .* repmat(t,2,1);
			t = (0:n2-1)/n2;
			BC = repmat(B,1,n2) + repmat(C-B,1,n2) .* repmat(t,2,1);
			t = (0:n3-1)/n3;
			CD = repmat(C,1,n3) + repmat(D-C,1,n3) .* repmat(t,2,1);
			t = (0:n4-1)/n4;
			DA = repmat(D,1,n4) + repmat(A-D,1,n4) .* repmat(t,2,1);
			
			points0 = [AB BC CD DA] ;
			
			theta = (0:nbPoints-1)/nbPoints * 2*pi;
			
			if dspl >= 1
				t0 = floor(n4/2);
				points = circshift(points0,[0,t0]);
				
				[points, tvec, avec, normal] = shape.C2boundary.rescale(points, theta, nbPoints, [], dspl);
			else
				tvec = [repmat((B-A),1,n1)/t1 repmat((C-B),1,n2)/t2 repmat((D-C),1,n3)/t3 repmat((A-D),1,n4)/t4]/2/pi ; % velocity vector
				
				rotation = [[0 1];[-1 0]] ;
				normal = rotation*tvec ;
				normal = normal./repmat(sqrt(normal(1,:).^2+normal(2,:).^2),2,1) ; % normal vector
				
				avec = zeros(2,nbPoints) ; % acceleration vector
				
				% shift the starting point to the mid-point of DA
				t0 = floor(n4/2);
				points = circshift(points0,[0,t0]);
				tvec = circshift(tvec,[0,t0]);
				normal = circshift(normal,[0,t0]);
			end
			
			if a==b
				name_str = 'Square';
			else
				name_str = 'Rectangle';
			end
			
			obj = obj@shape.C2boundary(points, tvec, avec, normal, [0,0]', name_str);
			
			obj.width = b;
			obj.height = a;
		end
		
		function obj = mtimes(obj, s)
			% Overload of the operator *
			obj = mtimes@shape.C2boundary(obj, s);
			obj.width = obj.width * s;
			obj.height = obj.height * s;
		end
		
	end
end
