classdef C2boundary
	% Abstract class for C2 smooth closed boundary.
	
	properties(SetAccess = protected)
		%% Quantites which are manually set
		
		points % coordinates of boundary points, an array of dimension 2 X nbPoints
		tvec % tangent vector
		avec % acceleration vector
		normal % outward normal vector
		center_of_mass % center of mass (run get_center_of_mass function only once)
		nbPoints % number of discrete boundary points (for runtime compatibility with the
		% functions like get.theta)
		
		name_str % name of the class
		
		%% Quantites which are automatically set
		
		theta % non tied-off parameterization between [0, 2pi)
		cpoints % complexification of the boundary points: points(1,:)+1i*points(2,:)
		diameter % (upper bound of) diameter of the shape, calculated from the center of mass
		tvec_norm % norm of the tangent vector
		pdirection % principle direction of the shape
		sigma % element of curve integration.
		% \int f(x) ds(x) = \int_0^(2pi) f(x(t)) |x'(t)| dt
		% ~ \sum_{n=1..nbPoints} f(points(n)) * sqrt(tvec_norm_square(n)) * 2pi/nbPoints
		% = \sum_{n=1..nbPoints} f(points(n)) * sigma(n)

        box % a minimal rectangular box [width, height] containing the shape (width: size in x-axis, height: size in y-axis)
	end
	
	methods
		function obj = C2boundary(points, tvec, avec, normal, com, nstr)
			% Initialization of C2boundary object:
			% Inputs:
			% points: coordinates of boundary point, a 2 X nbPoints array
			% tvec: tangent vectors of boundary point, a 2 X nbPoints array
			% avec: acceleration vectors, a 2 X nbPoints array
			% normal: outward normal vectors
			% com: center of mass of the boundary, optional (can be [])
			% nstr: name of the boundary, a string, optional
			
			obj.points = points;
			obj.tvec = tvec;
			obj.avec = avec;
			obj.normal = normal;
			
            % Check the curve
            flag = shape.C2boundary.check_sampling(points);
            if ~flag
                warning('Singularities found in the curve! Decrease the sampling step!');
            end
            
			if nargin > 4 && ~isempty(com)
				obj.center_of_mass = com;
			else
				obj.center_of_mass = shape.C2boundary.get_com(points, tvec, normal);
			end
			
			if nargin > 5
				obj.name_str = nstr;
			else
				obj.name_str = '';
			end
		end
		
		function val = get.box(obj)
			dd = obj.points - repmat(obj.center_of_mass, 1, obj.nbPoints);
			w = max(dd(1,:)) - min(dd(1,:));
			h = max(dd(2,:)) - min(dd(2,:));
			val = [w, h];
		end
		
		function val = get.pdirection(obj)
			D= obj.points - repmat(obj.center_of_mass, 1, obj.nbPoints);
			M = D*D';
			[P,~] = svd(M);
			X = P(:,1);
			t = atan(X(2)/X(1));
			val = [cos(t); sin(t)];
		end
		
		function val = get.nbPoints(obj)
			val = length(obj.points);
		end
		
		function val = get.diameter(obj)
			dd = obj.points - repmat(obj.center_of_mass, 1, obj.nbPoints);
			val = 2*max(sqrt(dd(1,:).^2 + dd(2,:).^2));
		end
		
		function val = get.theta(obj)
			val = 2*pi*(0:obj.nbPoints-1)/obj.nbPoints; % non tied-off
			
			% theta = 0:(2*pi/(nbPoints-1)):2*pi; % tied-off
		end
		
		function val = get.cpoints(obj)
			val = obj.points(1,:)+1i*obj.points(2,:);
		end
		
		function val = get.tvec_norm(obj)
			val = sqrt(obj.tvec(1,:).^2 + obj.tvec(2,:).^2);
		end
		
		function val = get.sigma(obj)
			val = 2*pi / obj.nbPoints * obj.tvec_norm;
		end
		
		%% Overloading of usual operators
		function obj = plus(obj, z0)
			% Overload of the operator +
			if ~isa(z0, 'double')
				error('Type error: only double value can be used for translation.');
			end
			if size(z0,1) ~= 2 || size(z0,2) ~= 1
				error('Size error: translation vector must have size (2,1)');
			end
			obj.points = obj.points + repmat(z0,1,obj.nbPoints);
			obj.center_of_mass = obj.center_of_mass + z0;
		end
		
		function obj = minus(obj, z0)
			% Overload of the operator -
			obj = plus(obj,-z0);
		end
		
		function obj = mtimes(obj, s)
			% Overload of the operator *
			if ~isa(s, 'double') || length(s) ~= 1 || s <= 0
				error('Type error: only positive double scalar can be used for scaling.');
			end
			
			obj.points = obj.points * s;
			obj.center_of_mass = obj.center_of_mass * s;
			obj.tvec = obj.tvec * s;
			obj.avec = obj.avec * s;
		end
		
		function obj = lt(obj, phi)
			% Redefine the < operator as the rotation of the boundary
			if ~isa(phi, 'double') || length(phi) ~= 1
				error('Type error: only double scalar can be used for rotation.');
			end
			
			rot = [[cos(phi), -sin(phi)]; [sin(phi), cos(phi)]]; % Rotation matrix
			
			obj.points = rot * obj.points;
			obj.center_of_mass = rot * obj.center_of_mass;
			obj.tvec = rot * obj.tvec;
			obj.normal = rot * obj.normal;
			obj.avec = rot * obj.avec;
		end
		
		%% Utility functions
		function obj1 = subset(obj, idx)
			% Return a subset of C2boundary object.
			if length(idx)>obj.nbPoints
				error('Value Error: the index of the subset is wrong!');
			end
			obj1 = obj; obj1.nbPoints = length(idx);
			obj1.points = obj.points(:, idx);
			obj1.tvec = obj.tvec(:, idx);
			obj1.avec = obj.avec(:, idx);
			obj1.normal = obj.normal(:, idx);
		end
		
		function val = isinside(obj, x)
			% check if a point is inside the ball B(center_of_mass, diameter/2)
			val = 1;
			if norm(x-obj.center_of_mass) >= obj.diameter/2
				val = 0;
			end
		end
		
		function val = isdisjoint(obj, E)
			% check if two boundaries are disjoint in the sense that the balls B(com1, radius1) and
			% B(com2, radius2) are disjoint.
			val = 1;
			d = norm(obj.center_of_mass - E.center_of_mass);
			if d <= (obj.diameter+E.diameter)/2
				val = 0;
			end
		end
		
		function plot(obj,varargin)
			plot(obj.points(1,:), obj.points(2,:),varargin{:}); hold on;
		end
		
		function [Sx, Sy, mask] = interior_mesh(obj, width, N)
			% Generate a square mesh centered at obj.center_of_mass and compute the mask of
			% the mesh points included in the interior of the domain.
			%
			
			N = 2*floor(N/2);
			z0 = obj.center_of_mass;
			[Sx, Sy, mask] = obj.boundary_off_mask(z0, width, N, width/N*2);
		end
		
		function [Sx, Sy, mask] = boundary_off_mask(obj, z0, width, N, epsilon)
			% Compute a binary mask of size N X N inside the square region (z0, width) so that the
			% mesh points on the boundary are turn off (0).
			% Inputs:
			% z0: center of the mesh
			% width: width of the mesh
			% N: number of points by side
			% Outputs:
			% Sx, Sy: mesh coordinates of boundary points
			% mask: binary mask
			mask = ones(N,N);
			
			sx = linspace(z0(1)-width/2, z0(1)+width/2, N);
			sy = linspace(z0(2)-width/2, z0(2)+width/2, N);
			
			[Sx, Sy] = meshgrid(sx, sy);
			Z = [Sx(:) Sy(:)]';
			
			for n=1:size(Z,2)
				dd = tools.dist_p2D(Z(:,n), obj.points);
				if dd<epsilon
					mask(n)=0;
				end
			end
			
			% % Another way: make the binary mask of the boundary
			% for n=1:obj.nbPoints
			%     x=obj.points(:,n)-z0;
			%     cidx = ceil(x(1)/dh) + N/2; % value from 1 to N
			%     ridx = N/2 - ceil(x(2)/dh);
			%     if (ridx > 0 && ridx <= N) && (cidx > 0 && cidx <= N)
			%         mask(ridx, cidx) = 1;
			%     end
			% end
			% % fill the inside of the mask
			% mask = imfill(mask, 4, 'holes'); % use matlab's imfill function
        end

        function obj1 = smooth(obj, epsilon, pos, width)
            % Smooth a piece of the boundary by convolution using a constant window of width
            % proportional to epsilon. The boundary to be smoothed is on [pos-width/2, pos+width/2]
        
            epsilon = mod(epsilon, 1);
            hwidth = ceil(epsilon * obj.nbPoints);
            
            if nargin < 3 
                p1 = tools.convfix(obj.points(1,:), hwidth);
                p2 = tools.convfix(obj.points(2,:), hwidth);
            else
                pos = mod(pos, 1);
                width = mod(width, 1);
                idx = floor(pos * obj.nbPoints);
                Lt = max(1, floor(obj.nbPoints * width)/2);

                q1 = tools.convfix(obj.points(1,idx-Lt:idx+Lt), hwidth);
                q2 = tools.convfix(obj.points(2,idx-Lt:idx+Lt), hwidth);
                p1 = [obj.points(1, 1:idx-Lt-1), q1, obj.points(1,idx+Lt+1:end)];
                p2 = [obj.points(2, 1:idx-Lt-1), q2, obj.points(2,idx+Lt+1:end)];
            end
            D = [p1; p2]; 
            N=length(p1); theta=2*pi*(0:N-1)/N;
            [D1, tvec1, avec1, normal1] = shape.C2boundary.rescale(D, theta, obj.nbPoints, obj.box);
            obj1 = shape.C2boundary(D1, tvec1, avec1, normal1, [], obj.name_str);
        end
        
		function obj1 = dammage(obj, epsilon, pos, width)
			% Apply a dammage to the boundary by linking boundary points.
        end
        
		function obj1 = global_perturbation(obj, epsilon, p, n)
			% Perturb and smooth globally a boundary:
			% D_epsilon(t) = D(t)(1+epsilon*cos(p*t)*normal(t))
			%
			% Inputs:
			% epsilon: strength of the perturbation
			% p: periodicity of the perturbation
			% n: strength of the smooth filter, integer
			%
			% Output: a new C2boundary object
			
			D = obj.points + epsilon*repmat(cos(p*obj.theta),2,1) .* obj.normal;
			
			% Smooth the bounary
			if n>0
				k = ones(1,n)/n;
				mode = 'same';
				
				M = length(D);
				d1 = [D(1,:), D(1,:), D(1,:)];
				d2 = [D(2,:), D(2,:), D(2,:)];
				
				d1 = conv(d1, k, mode);
				d2 = conv(d2, k, mode);
				D = [d1(M+1:2*M); d2(M+1:2*M)];
			end
			
			[D1, tvec1, avec1, normal1] = shape.C2boundary.rescale(D, obj.theta, obj.nbPoints, obj.box);
			obj1 = shape.C2boundary(D1, tvec1, avec1, normal1, [], obj.name_str);
			%            obj1 = obj.recreat(D, obj.theta, obj.nbPoints);
		end
		
		function obj1 = local_perturbation(obj, epsilon, pos, width, theta)
			% Perturb locally a boundary:
			% D_epsilon(t) = D(t)(1+epsilon*h(t)*rot(theta)*normal(t))
			% h(t) is a compactly supported function, rot(theta) is the
			% rotation matrix with angle theta
			%
			% Inputs:
			% epsilon: strength of the perturbation
			% pos: 0<=pos<=1 is the center of the support of h
			% width: 0<=width<=1 is the width of the support of h, and
			% [pos-width/2, pos+width/2] is the support
			% theta: an angle, optional
			%
			% Output: a new C2boundary object
			
			if width > 0.5
				error('Wrong value of width, must be smaller than 0.5!');
			end
			
			pos = mod(pos, 1);
			width = mod(width, 1);
			idx = floor(pos * obj.nbPoints);
			
			% h = @(t)(exp(-1./(1-t.^2))); % a C^infty function with compact support
			h = @(t)(exp(-10*t.^2));
			Lt = max(1, floor(obj.nbPoints * width));
			toto = [h(linspace(-1,1,Lt)), zeros(1,obj.nbPoints-Lt)];
			H = circshift(toto, [0, idx-floor(Lt/2)]);
						
			if nargin < 5 || isempty(theta)
				rot = eye(2);
			else
				rot = [cos(theta), sin(-theta); sin(theta), cos(theta)];
			end
			D = obj.points + epsilon * repmat(H, 2, 1) .* (rot * obj.normal);
			% figure; plot(D(1,:), D(2,:));
			
			[D1, tvec1, avec1, normal1] = shape.C2boundary.rescale(D, obj.theta, obj.nbPoints);
			obj1 = shape.C2boundary(D1, tvec1, avec1, normal1, [], obj.name_str);
		end
	end
	
	methods(Access = protected)
		function val = get_center_of_mass(obj)
			val = shape.C2boundary.get_com(obj.points, obj.tvec, obj.normal);
		end
		
	end
	
	methods(Static)
		[tvec,avec,normal] = boundary_vec(D, t)
		
		function [D, tvec, avec, normal] = rescale(D0, theta0, nbPoints, nsize)
			% Compute all variables related to the boundary from the boundary points D0 and the
			% parameterization theta0. The new boundary will be reinterpolated with nbPoints, and
			% optionally rescaled to fit the rectangle of size nsize=[width, height].
			
			if nargin == 4 && ~isempty(nsize) 
				minx = min(D0(1,:)) ;
				maxx = max(D0(1,:)) ;
				miny = min(D0(2,:)) ;
				maxy = max(D0(2,:)) ;
				
				z0 = [(minx+maxx)/2; (miny+maxy)/2];
				D0 = [(D0(1,:)-z0(1))*(nsize(1)/(maxx-minx)); (D0(2,:)-z0(2))*(nsize(2)/(maxy-miny))];
			end
			
			theta = (0:nbPoints-1)/nbPoints*2*pi;
			D = interp1(theta0, D0', theta, 'spline');
			D = D';
			
			% high order bounary informations
			[tvec,avec,normal] = shape.C2boundary.boundary_vec(D, theta);
			
			% obj.points = D;
			% obj.tvec = tvec;
			% obj.avec = avec;
			% obj.normal = normal;
			% obj.name_str = name_str;
			% obj.center_of_mass = obj.get_center_of_mass(); % run only once, for speed
		end
		
		%        D = rescale_shape(D0, w, h, nbPoints)
		
		function val = get_com(points, tvec, normal)
			% Calculate the center of mass of a shape by the Stokes formula
			% INPUTS:
			% points: coordinates of boudary points
			% normal: normal vector of the boudary
			% sigma: boundary integral elements
			
			nbPoints = length(points);
			tvec_norm = sqrt(tvec(1,:).^2 + tvec(2,:).^2);
			sigma = 2*pi / nbPoints * tvec_norm;
			mass = (sum(points(1,:).*normal(1,:).*sigma) + sum(points(2,:).*normal(2,:).*sigma))/2;
			
			Cx = sum(1/2 * (points(1,:).^2) .* normal(1,:).*sigma);
			Cy = sum(1/2 * (points(2,:).^2) .* normal(2,:).*sigma);
			
			val = [Cx,Cy]' / mass;
			%C = (Cx+j*Cy) / mass;
        end
        
        function val = check_sampling(points)
            % A C^1 parameterized simple curve f(t) (t is the parameter) must satisfy (by Taylor expansion)
            % <f(t_n+1)-f(t_n), f(t_n)-f(t_n-1)> > 0, for sufficiently
            % small sampling step dt = t_n - t_n-1.
            %
            % This function check this condition.
            
            nbPoints = size(points, 2);

            val = 1;

            for p=1:nbPoints
                x = points(:,p);
                
                if p==1
                    y=points(:,nbPoints);
                    z=points(:,p+1);
                elseif p==nbPoints
                    y=points(:,p-1);
                    z=points(:,1);
                else
                    y=points(:,p-1);
                    z=points(:,p+1);                    
                end
                
                toto = (z-x)'*(x-y)/(norm(y-x)*norm(z-x));
                if toto <= 0
                    val = 0;
                end
            end
        end

        function val = check_C1simplecurve(points)
            % A C^1 curve is simple means it is a) connected, b) does not cross
            % itself and c) ends at the starting point. This function check
            % the condition b).
            
           % TODO
           val = 0;
        end
        
	end
end

