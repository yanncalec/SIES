classdef Banana < shape.C2boundary
	% Class for banana-shaped object, used mainly in ElectricFish class
	
	properties(SetAccess = protected)
		axis_a = 1; % length of the semi-major axis
		axis_b = 1; % length of the semi-minor axis
		center;
		center_curvature;
	end
	
	methods
		function obj = Banana(a,b,center,xc,yc,nbPoints)
			% This function creates a curved ellipse, which looks like a banana
			% INPUT : - a,b = semi axis of the ellipse
			%         - center = coordinates of the ellipse's center
			%         - xc,yc = coordinates of the center of curvature
			%         - N = number of points
			
			if a<b
				error(['Value error: the semi-major axis must be longer than the ' ...
					'semi-minor one.']);
			end
			
			x0 = center(1); y0 = center(2);
			
			% La ligne de courbure de l'ellipse
			R = sqrt((xc-x0)^2 + (yc-y0)^2) ;
			theta0 = atan2(y0-yc,x0-xc) ;
			alpha = a/R ;
			
			theta = 2*pi*(0:nbPoints-1)/nbPoints;
			t = theta0+alpha*cos(theta) ;
			points = [xc + (R+ b*sin(theta)).*cos(t) ; yc + (R+ b*sin(theta)).*sin(t)]  ;
			
			tvec = [ b*cos(theta).*cos(t) + alpha*sin(theta).*(R+b*sin(theta)).*sin(t) ; ...
				b*cos(theta).*sin(t) - alpha*sin(theta).*(R+b*sin(theta)).*cos(t) ] ;
			
			rotation = [[0 -1];[1 0]] ;
			normvec = rotation*tvec ;
			normal = normvec./repmat(sqrt(normvec(1,:).^2+normvec(2,:).^2),2,1) ;
			
			avec = [-b*sin(theta).*cos(t) + alpha*b*cos(theta).*sin(t) + alpha*cos(theta).*(R+b*sin(theta)).*sin(t) + ...
				alpha*b*cos(theta).^2.*sin(t)- alpha^2*sin(theta).^2.*(R+b*sin(theta)).*cos(t) ; ...
				-b*sin(theta).*cos(t) - alpha*b*sin(theta).*cos(t) - alpha*cos(theta).*(R+b*sin(theta)).*cos(t) + ...
				-alpha*b*cos(theta).^2.*cos(t)- alpha^2*sin(theta).^2.*(R+b*sin(theta)).*sin(t)] ;
			
			obj = obj@shape.C2boundary(points, tvec, avec, normal, [], 'Banana');
			
			obj.center = center(:);
			obj.center_curvature = [xc; yc];
			obj.axis_a = a;
			obj.axis_b = b;
		end
		
		%% Overloading of usual operators
		function obj = plus(obj, z0)
			obj = plus@shape.C2boundary(obj, z0);
			obj.center = z0+obj.center;
			obj.center_curvature = z0+obj.center_curvature;
		end
		
		function obj = mtimes(obj, s)
			obj = mtimes@shape.C2boundary(obj,s);
			obj.center = s*obj.center;
			obj.center_curvature = s*obj.center_curvature;
			obj.axis_a = obj.axis_a * s;
			obj.axis_b = obj.axis_b * s;
		end
		
		function obj = lt(obj, phi)
			obj = lt@shape.C2boundary(obj,phi);
			rot = [[cos(phi), -sin(phi)]; [sin(phi), cos(phi)]]; % Rotation matrix
			obj.center = rot*obj.center;
			obj.center_curvature = rot*obj.center_curvature;
		end
	end
end