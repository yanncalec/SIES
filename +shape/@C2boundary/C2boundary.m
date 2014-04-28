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
    end
    
    methods        
        function obj = C2boundary(points, tvec, avec, normal, com, nstr)
            obj.points = points;
            obj.tvec = tvec;
            obj.avec = avec;
            obj.normal = normal;

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
    
        function obj1 = global_perturbation(obj, epsilon, p, n)
        % Perturb and smooth globally a boundary:
        % D_epsilon(t) = D(t)(1+epsilon*cos(t)*normal(t))
        %
        % Inputs:
        % p: periodicity of the pertubation
        % n: strength of the smooth filter, integer
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

            [D1, tvec1, avec1, normal1] = shape.C2boundary.rescale(D, obj.theta, obj.nbPoints);
            obj1 = shape.C2boundary(D1, tvec1, avec1, normal1, [], obj.name_str);
            %            obj1 = obj.recreat(D, obj.theta, obj.nbPoints);
        end
        
        % function obj1 = local_perturbation(obj)
        %     error('Not implemented');
        % end
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
            
            if ~isempty(nsize)
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
    end
end

% properties
%     cnd = 1; % conductivity of the shape, a positive number.
%              % Set cnd small for conductive object, and cnd large for resistive object
%     pmtt = 1; % permittivity of the shape
%     pmeb = 1; % mu, permeability of the shape                
% end

% %% Get and Set functions

% function obj = set.cnd(obj, val)
%     if ~isreal(val) || val == 1
%         error('Value error: conductivity must be a real number different to 1');
%     else
%         obj.cnd = val;
%     end
% end

% function obj = set.pmtt(obj, val)
%     if ~isreal(val) || val < 0
%         error('Value error: permittivity must be a positive real number');
%     else
%         obj.pmtt = val;
%     end
% end

% function obj = set.pmeb(obj, val)
%     if ~isreal(val) || val < 0
%         error('Value error: permeability must be a positive real number');
%     else
%         obj.pmeb = val;
%     end
% end


%         function val = subsref(obj, idx)
%             val = builtin('subsref', obj.points, idx);
%         end

%         function val = get.nbPoints(obj)
%             val = size(obj.points,2);
%         end


%         function val = get_center_of_mass(obj)
%             % Calculate the center of mass of a shape by the Stokes formula
%             mass = (sum(obj.points(1,:).*obj.normal(1,:).*obj.sigma) + sum(obj.points(2,:).*obj.normal(2,:).*obj.sigma))/2;
%             
%             Cx = sum(1/2 * (obj.points(1,:).^2) .* obj.normal(1,:).*obj.sigma);
%             Cy = sum(1/2 * (obj.points(2,:).^2) .* obj.normal(2,:).*obj.sigma);
%             
%             val = [Cx,Cy]' / mass;
%             %C = (Cx+j*Cy) / mass;
%         end


% function obj = t_s_r(obj, z0, s, phi)
% % Apply (by order) the rotaion, scaling and translation on boundary and update its related boundary
% % parameterizations.
% % Inputs:
% % z0: translation vector
% % phi: angle of rotaion 
% % s: scaling constant, s>0

%     obj = ((obj < phi) * s) + z0;
% end

% function obj = r_s_t(obj, phi, s, z0)
% % Apply (by order) the translation, scaling and rotaion on a boundary and update its related
% % boundary parameterizations.            
% % Inputs:
% % phi: angle of rotaion
% % s: scaling constant, s>0
% % z0: translation vector

%     obj = rotation((obj + z0) * s, phi);
% end

% function val = kvalH(obj, freq)
% % the frequency dependant wave number k for Helmholtz eq
%     val = freq*sqrt(obj.pmtt*obj.pmeb);
% end

% function val = wavenb(obj, freq)
% % the frequency dependant wave number k for Helmholtz eq
%     val = freq*sqrt(obj.pmtt*obj.pmeb);
% end

% function val = kappa(obj, freq)
% % the frequency dependant complex conductivity
% % kappa = conductivity + 1i * frequency * permittivity
%     val = obj.cnd + 1i*obj.pmtt * freq;
% end

% function val = lambda(obj, freq)
% % the frequency dependant complex contrast
% % lambda = (kappa+1)/2/(kappa-1)
%     kappa = obj.kappa(freq);
%     val = (kappa+1)./(kappa-1)/2 ;
% end

% function obj1 = dwspl(obj, df)
% % Downsample the boundary by a factor df and return a
% % C2boundary object.
%     if mod(obj.nbPoints, df) || df >= obj.nbPoints || df <= 0
%         error('Value Error: down-sampling factor is wrong!');
%     end
%     obj1 = obj; obj1.nbPoints = obj.nbPoints / df;
%     obj1.points = obj.points(:, 1:df:end);
%     obj1.tvec = obj.tvec(:, 1:df:end);
%     obj1.avec = obj.avec(:, 1:df:end);
%     obj1.normal = obj.normal(:, 1:df:end);
% end

% function obj1 = subset(obj, idx)
% % Return a subset of C2boundary object.
%     if length(idx)>obj.nbPoints
%         error('Value Error: the index of the subset is wrong!');
%     end
%     obj1 = obj; obj1.nbPoints = length(idx);
%     obj1.points = obj.points(:, idx);
%     obj1.tvec = obj.tvec(:, idx);
%     obj1.avec = obj.avec(:, idx);
%     obj1.normal = obj.normal(:, idx);
% end




% %% Generalized Polarization Tensors related
% function [Mc, M] = GPT(D, ord, freq)
% % Calculate the GPT of a shape by numerical integration
%     if nargin<3
%         freq = 0;
%     end
%     [Mc, M] = GPT.theoretical_GPT(D.points, D.tvec, D.normal, D.avec, D.sigma, D.kappa(freq), ord);            
% end

% function M = CGPT(D, ord, freq)
% % Calculate the CGPT of a shape by numerical integration
%     if nargin<3
%         freq = 0;
%     end
%     X = GPT.theoretical_CGPT(D.points, D.tvec, D.normal, D.avec, D.sigma, D.kappa(freq), ord);
%     M = GPT.CGPT(X);
% end

% function W = SCT(D, ord, freq, pmtt_bg, pmeb_bg)
% % Calculate the SCT (scattering coefficients) of a shape by numerical integration
% % INPUTS:
% % ord: order of SCT matrix. The SCT matrix has dimension (2*ord+1)^2
% % freq: working frequencies
% % pmtt_bg: permettivity of the background
% % pmeb_bg: permeability of the background
% % OUPTUP:
% % W: a 3d array of dimension (2*ord+1) X (2*ord+1) X length(freq)

%     Nf = length(freq);
%     W = zeros(2*ord+1, 2*ord+1, Nf);

%     for n=1:Nf
%         W(:,:,n) = tools.Helmholtz.theoretical_SCT(D, ord, freq(n), pmtt_bg, pmeb_bg);
%     end                
%     W = squeeze(W);
% end

% %% WPT
% function M = WPT_std(D, WF, freq, nbScl)
% % Calculate the Wavelet Polarization Tensor of a shape using
% % the standard representation. 
% % Inputs:
% % WF: an object of the directional wavelet frame
% % freq: frequency
% % nbScl: number of scales

%     if nargin<3
%         freq = 0;
%     end

%     if ~isa(WF, 'wavelet.Frame')
%         error('An object of wavelet.Frame is needed.');
%     end

%     M = WF.WPT_std(D, freq, nbScl);
% end


% function obj1 = rescale(obj, w, h, nbPoints)
% % Fit the boundary to a rectangle of width w and height h, and
% % re-interpolate it using nbPoints. The output is the new boundary.

%     minx = min(obj.points(1,:)) ;
%     maxx = max(obj.points(1,:)) ;
%     miny = min(obj.points(2,:)) ;
%     maxy = max(obj.points(2,:)) ;

%     z0=[(minx+maxx)/2; (miny+maxy)/2];

%     D0 = [(obj.points(1,:)-z0(1))*(w/(maxx-minx));...
%           (obj.points(2,:)-z0(2))*(h/(maxy-miny))] ; 

%     t0 = (0:nbPoints-1)/nbPoints*2*pi;

%     obj1 = shape.C2boundary(D0, t0, obj.name_str);
% end

% function obj1 = reinterpolate(obj, nbPoints)
% % Reinterpolate a boundary with nbPoints

%     t = (0:nbPoints-1)/nbPoints*2*pi;

%     D = interp1(obj.theta, obj.points', t, 'spline') ;
%     D = D';

%     obj1 = C2boundary(D, obj.name_str);
% end

        % function obj1 = recreat(obj, D0, theta0, nbPoints, nsize)
        % % Recreat a C2boundary object from the current object using new boundary points D0 and
        % % parameterization theta0. The new boundary will be reinterpolated with nbPoints, and
        % % optionally rescaled to fit the rectangle of size nsize=[width, height].
        
        %     obj1 = obj; 
        %     [D, tvec, avec, normal] = shape.C2boundary.rescale(D0, theta0, nbPoints, nsize);
        
        %     obj1.points = D;
        %     obj1.nbPoints = nbPoints;
        %     obj1.tvec = tvec;
        %     obj1.avec = avec;
        %     obj1.normal = normal;
        %     obj1.center_of_mass = obj1.get_center_of_mass();
        % end            
        
