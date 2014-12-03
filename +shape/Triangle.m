classdef Triangle < shape.C2boundary
	% Class for isosceles triangle
	
	properties
		lside % length of the equal side
		angl % angle between the two equal sides
	end
	
	methods
		function obj = Triangle(a, angl, nbPoints, dspl)
			% INPUT :
			% a: length of the equal side
			% angl: angle between the two equal sides
			% nbPoints: number of discretization points
			% dspl: down-sampling factor for smoothing the corners
			
			if nargin<4
				dspl = 10;
			end
			
			h = a*cos(angl/2); % height of the triangle
			b = a*sin(angl/2); % 2b is the length of the 3rd side
			
			t1 = a/(a+b)/2; t2 = b/(a+b); t3 = t1;
			n1 = floor(t1*nbPoints); n2 = floor(t2*nbPoints); n3 = nbPoints-n1-n2;
			
			A = [ 0 ; 2/3*h ] ; B = [ - b ; -h/3 ] ; C = [ b ; -h/3 ] ;
			
			t = (0:n1-1)/n1;
			AB = repmat(A,1,n1) + repmat(B-A,1,n1) .* repmat(t,2,1);
			t = (0:n2-1)/n2;
			BC = repmat(B,1,n2) + repmat(C-B,1,n2) .* repmat(t,2,1);
			t = (0:n3-1)/n3;
			CA = repmat(C,1,n3) + repmat(A-C,1,n3) .* repmat(t,2,1);
			% CA = (1-t)*repmat(C,1,M) + t*repmat(A,1,M);
			
			points0 = [AB BC CA] ;
			theta = (0:nbPoints-1)/nbPoints * 2*pi;
			
			if dspl >= 1
				t0 = floor(n3/2);
				points = circshift(points0,[0,t0]);
				
				[points, tvec, avec, normal] = shape.C2boundary.rescale(points, theta, nbPoints, [], dspl);
			else
				tvec = [repmat((B-A),1,n1)/t1 repmat((C-B),1,n2)/t2 repmat((A-C),1,n3)/t3]/2/pi ; % velocity vector
				
				rotation = [[0 1];[-1 0]] ;
				normal = rotation*tvec ;
				
				normal = normal./repmat(sqrt(normal(1,:).^2+normal(2,:).^2),2,1) ;
				avec = zeros(2,nbPoints) ; % acceleration vector
				
				% shift the starting point to the mid-point of CA
				t0 = floor(n3/2);
				points = circshift(points0,[0,t0]);
				tvec = circshift(tvec,[0,t0]);
				normal = circshift(normal,[0,t0]);
				
			end
			
			obj = obj@shape.C2boundary(points, tvec, avec, normal, [0,0]', 'Triangle'); % the triangle is centered (center of the mass) at the origine by construction
			% smooth out the singularity
			obj.lside = a;
			obj.angl = angl;
		end
		
		function obj = mtimes(obj, s)
			% Overload of the operator *
			obj = mtimes@shape.C2boundary(obj, s);
			obj.lside = obj.lside * s;
		end
		
	end
end

