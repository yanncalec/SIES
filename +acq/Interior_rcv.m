classdef Interior_rcv < acq.mconfig
% Acquisition system with sources uniformly distributed on a circle, and receivers uniformly on a
% Cartesian grid in the interior of the circular ROI.
    
    properties(SetAccess = protected)
        radius % radius of the source circle

        ROI % Square region of interest, wavelet.ROI class
        dim % [row, col] dimension of the cartesian grid
        
        dh % sampling step on the grid
        rangex % a vector indicating interger index range of the sampling points in x-axis of ROI
        rangey % a vector indicating interger index range of the sampling points in y-axis of ROI
        meshpoints
        mask % binary mask for active receivers
        regdist=1e-2 % minimal distance between the singularity of Green function and the boundary
    end
    
    methods
        function obj = Interior_rcv(Z, width, N, Rs, Ns, modes, D)
        % Inputs:
        % Z: center of the measurement circle
        % width: width of the ROI
        % N: number of sampling points by side of the ROI
        % Rs: radius of the source circle
        % Ns: number of sources per arc 
        % modes: viewmode, see acq.Concentric
        % D: a boundary (array of coordinates), optional. If given, all receivers close to the boundary will be removed.

            if nargin < 7
                D=[];
            end
            
            if nargin < 6
                modes = [1, 2*pi, 2*pi];
            end
            
            obj.center = Z;
            obj.radius = Rs;
            obj.dh = width/N;
            obj.ROI = wavelet.ROI(obj.center, width, width);            
            
            Ngs = modes(1); thetas = modes(2); aovs = modes(3);
            
            toto = acq.src_rcv_circle(Ngs, Ns, Rs, Z, thetas, aovs);
            obj.src_prv = {toto};
            
            % obj.aov_src = min(2*pi, max(aovs, Ngs*thetas));
            
            tx0 = ceil((Z(1)-width/2)/obj.dh);
            tx1 = min(N+tx0-1, floor((Z(1)+width/2)/obj.dh));
            ty0 = ceil((Z(2)-width/2)/obj.dh);
            ty1 = min(N+ty0-1, floor((Z(2)+width/2)/obj.dh));
            obj.rangex = [tx0, tx1];
            obj.rangey = [ty0, ty1];

            obj.dim = [obj.rangey(2)-obj.rangey(1)+1, obj.rangex(2)-obj.rangex(1)+1];

            obj.mask = ones(obj.dim);
            
            if ~isempty(D)
                Nx = obj.dh*(obj.rangex(1):obj.rangex(2));
                Ny = obj.dh*(obj.rangey(1):obj.rangey(2));
                [Sx, Sy] = meshgrid(Nx, Ny);
                rcv = [Sx(:) Sy(:)]';
                
                for n=1:prod(obj.dim) 
                    dd = tools.dist_p2D(rcv(:,n), D);
                    if dd<obj.regdist || norm(rcv(:,n)-obj.center) > obj.radius
                       obj.mask(n)=0;
                    end
                end
            end

            toto = obj.meshpoints;
            obj.rcv_prv = {toto};
            
            obj.Ng = 1; % Full view mode only
        end
        
        function val = get.meshpoints(obj)
            Nx = obj.dh*(obj.rangex(1):obj.rangex(2));
            Ny = obj.dh*(obj.rangey(1):obj.rangey(2));
            [Sx, Sy] = meshgrid(Nx, Ny);
            Sx = Sx(find(obj.mask));
            Sy = Sy(find(obj.mask));
            
            val = [Sx(:) Sy(:)]';
        end

        % function val = rcv(obj, n)
        %     if n>obj.Ns
        %         error('Source index out of range');
        %     end
        %     val = obj.rcv_prv{1};
        % end
        
        function plot(obj, varargin)            
            src = obj.src_prv{1};
            plot(src(1,:), src(2,:), 'x'); hold on;
            rcv = obj.rcv_prv{1};
            plot(rcv(1,:), rcv(2,:), '.');            
            plot(obj.center(1), obj.center(2), 'r*');
        end        
    end
end

