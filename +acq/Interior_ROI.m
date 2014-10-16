classdef Interior_ROI < acq.mconfig
	% Acquisition system with coincided sources/receivers distributed uniformly
	% on a Cartesian grid inside a squared ROI.
	
	properties(SetAccess = protected)
		ROI % object of wavelet.ROI class
		dim % [row, col] dimension of the cartesian grid
		
		dh % sampling step of the grid
		rangex % a vector indicating interger index range of the sampling points in x-axis of ROI
		rangey % a vector indicating interger index range of the sampling points in y-axis of ROI
		meshpoints
		mask % binary mask for active receivers
		
		regdist=1e-3 % minimal distance between the singularity of Green function and the boundary
	end
	
	methods
		function obj = Interior_ROI(Z, width, N, D)
			% INPUTS:
			% -Z: center of the ROI
			% -width: width of the ROI
			% -N: number of sampling points by side
			% -D: a reference boundary inside ROI, in order to remove the sources/receivers
			% located exactly on the boundary of D. Empty by default.
			
			if nargin < 4
				D=[];
			end
			
			obj.center = Z;
			obj.dh = width/N;
			obj.ROI = wavelet.ROI(obj.center, width, width);
			
			tx0 = ceil((Z(1)-width/2)/obj.dh);
			tx1 = min(N+tx0-1, floor((Z(1)+width/2)/obj.dh));
			ty0 = ceil((Z(2)-width/2)/obj.dh);
			ty1 = min(N+ty0-1, floor((Z(2)+width/2)/obj.dh));
			obj.rangex = [tx0, tx1];
			obj.rangey = [ty0, ty1];
			
			obj.dim = [obj.rangey(2)-obj.rangey(1)+1, obj.rangex(2)-obj.rangex(1)+1];
			
			% obj.fixed_rcv = 1;
			% obj.coincided_src_rcv = true;
			
			obj.mask = ones(obj.dim);
			
			if ~isempty(D)
				Nx = obj.dh*(obj.rangex(1):obj.rangex(2));
				Ny = obj.dh*(obj.rangey(1):obj.rangey(2));
				[Sx, Sy] = meshgrid(Nx, Ny);
				src = [Sx(:) Sy(:)]';
				
				for n=1:prod(obj.dim)
					dd = tools.dist_p2D(src(:,n), D);
					if dd<obj.regdist
						obj.mask(n)=0;
					end
				end
			end
			
			obj.src_prv = {obj.meshpoints};
			obj.rcv_prv = {obj.meshpoints};
			obj.Ng = 1;
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
			plot(src(1,:), src(2,:), '.'); hold on;
			plot(obj.center(1), obj.center(2), 'r*');
		end
	end
end

