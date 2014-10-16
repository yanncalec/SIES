classdef Fish_circle < acq.mconfig
	% Circular acquisition for electric fish.
	%
	% This configuration is useful for electric fish class |PDE.ElectricFish|.
	% The fish has one eletric organ (source) and many receivers on its body, and it
	% swims in a circular trajectory while taking measurments. The position of
	% receivers depends on that of the source since they move with the fish. The
	% dipole source changes its direction with the fish's position.
	
	properties(SetAccess = protected)
		dipole_prv % directions of dipole source, a cell of 2D vectors
	end
	
	properties(SetAccess = protected)
		Omega0 % a copy of the original fish body
		impd % impedance of the fish's skin
		dipole0 % direction the dipole source in Omega0, a 2D vector
		eorgan0 % position of the eletric organ in Omega0, a 2D vector
		
		aov % total angle of view effectively covered by sources
		radius % radius of measurement circle for sources
		
		angl % rotation angle of the fish's body at different source positions
		fpos % fish's position. Be careful to not confuse fpos with src, which is the position of the
		% electric organ.
		
		idxRcv % index indicating the presence of receiver
	end
	
	methods
		function obj = Fish_circle(Omega, idxRcv, Z, Rs, Ns, aov, eorgan0, dipole0, ...
				d0, impd)
			% obj = Fish_circle(Omega, idxRcv, Z, Rs, Ns, aov, eorgan0, dipole0, d0, impd)
			% Omega: body of the fish, a C2boundary object
			% idxRcv: index of active receivers. If idxRcv is empty, then all
			%         boundary points will be receivers.
			% Z: center of the measurement circle
			% Rs: radius of the measurement circle
			% Ns: number of sources
			% aov: angle of view covered by all sources
			% eorgan0: position of reference of the electric organ (optional)
			% dipole0: direction of reference of the dipole source (optional)
			% d0: offset of the dipole source wrt the center of mass (optional)
			% impd: impedance of the skin (optional)
			
			if ~isa(Omega, 'shape.C2boundary')
				error('Type Error: the domain must be an object of C2boundary');
			end
			obj.Omega0 = Omega;
			
			if nargin < 10 || isempty(impd)
				obj.impd = 0; % by default the impedance of skin is 0
			else
				obj.impd = impd;
			end
			
			if nargin < 9 || isempty(d0)
				d0 = [0,0]'; % no offset by default
			end
			
			% Dipole of reference
			if nargin < 8 || isempty(dipole0)
				% by default the dipole direction is the shape's principle direction
				obj.dipole0 = Omega.pdirection;
			else
				obj.dipole0 = dipole0(:)/norm(dipole0);
			end
			
			% Electric organ of reference
			if nargin < 7 || isempty(eorgan0)
				% by default the electric organ is at the center of the mass, with an offset d0
				obj.eorgan0 = Omega.center_of_mass + diag(d0) * obj.dipole0 * Omega.diameter/2;
			else
				obj.eorgan0 = eorgan0;
			end
			
			% Body of reference
			if isa(Omega, 'shape.Banana') % special treatment for the banana-shaped fish
				xc = Omega.center_curvature(1);
				yc = Omega.center_curvature(2);
				x0 = Omega.center(1); y0 = Omega.center(2);
				R = sqrt((xc-x0)^2 + (yc-y0)^2) ;
				% Thomas Boulier's version
				%                 theta0 = atan2(y0-yc,x0-xc) ;
				%                 a = Omega.axis_a; b = Omega.axis_b;
				%                 alpha = 2*a/R ;
				%                 angle_source = theta0+alpha/2*d0(1) ;
				%                 source_position = [xc  yc] + (R+b*d0(2))*[cos(angle_source)  sin(angle_source)]  ;
				%
				%                 theta_p = 0; % angle between the line of the fish and the dipole source
				%                 obj.dipole0 = [cos(angle_source-pi/2+theta_p) sin(angle_source-pi/2+theta_p)]' ;
				%                 obj.eorgan0 = source_position(:);
				
				alpha = 2*d0(1) * atan(max(Omega.axis_a, Omega.axis_a)/R);
				obj.dipole0 = [-sin(alpha); cos(alpha)];
				obj.eorgan0 = R*[cos(alpha); sin(alpha)];
			end
			
			obj.aov = aov;
			obj.center = Z(:);
			
			% Compute the positions of the fish's body
			if isa(Omega, 'shape.Banana') % special treatment for the banana-shaped fish
				obj.radius = 0; % The rotation of the babana fish is special. In this case the radius doesn't make
				% sense.
			else
				obj.radius = Rs;
			end
			
			[obj.fpos, obj.angl] = acq.src_rcv_circle(1, Ns, obj.radius, Z, aov, 2*pi);
			
			% Active receivers
			if isempty(idxRcv)
				obj.idxRcv = 1:Omega.nbPoints;
			else
				obj.idxRcv = idxRcv;
			end
			
			% Compute positions of receivers
			src_toto = zeros(2,Ns);
			obj.Ng = Ns; % number of the group
			
			for s=1:Ns
				theta = obj.angl(s);
				Body = (obj.Omega0 < theta) + obj.fpos(:,s);
				obj.rcv_prv{s} = Body.points(:, obj.idxRcv); % receiver positions corresponding to the s-th source
				
				% if s==1 % explicitely store the receivers' position only for the first source
				%     Body = (obj.Omega0 < theta) + obj.fpos(:,s);
				%     obj.rcv_prv{1} = Body.points(:, obj.idxRcv); % receiver positions corresponding to the s-th source
				% end
				
				rot = [cos(theta) -sin(theta); sin(theta) cos(theta)];
				obj.dipole_prv{s} = rot*obj.dipole0; % the direction of s-th dipole
				obj.src_prv{s} = obj.fpos(:,s) + rot*obj.eorgan0; % the position of s-th electric organ
			end
		end
		
		function val = all_dipole(obj)
			% Export all dipoles
			val = [];
			for g=1:obj.Ng
				val = [val obj.dipole_prv{g}];
			end
		end
		
		function val = dipole(obj, n)
			% n-th dipole direction
			if n>obj.Ns_total || n<1
				error('Source index out of range');
			end
			val = obj.dipole_prv{n};
		end
		
		function val = Bodies(obj, n)
			% Return a C2boundary object containning the fish body at the
			% n-th position. For a position of sources, we rotate and
			% translate the body of the fish.
			if n>obj.Ns_total || n<1
				error('Source index out of range.');
			end
			val = (obj.Omega0 < obj.angl(n)) + obj.fpos(:,n);
		end
		
		function plot(obj, idx, varargin)
			
			if nargin<2 || isempty(idx)
				idx = 1:obj.Ns_total;
			end
			
			for n=1:length(idx)
				s = idx(n);
				rcv = obj.rcv(s);
				plot(rcv(1,:), rcv(2,:), '.', varargin{:}); hold on;
				Omega = obj.Bodies(s); plot(Omega);
				
				src = obj.src(s);
				dp = obj.dipole(s)*Omega.diameter*0.25;
				plot(src(1), src(2), 'go', varargin{:});
				
				quiver(src(1), src(2), dp(1), dp(2), 'MarkerSize',.1, varargin{:}); % dipole vector
			end
			plot(obj.center(1), obj.center(2), 'r*', varargin{:});
		end
	end
	
end
