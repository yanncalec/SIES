classdef Small_Inclusions < handle
	% Abstract class for several PDEs involving small inclusions in a homogeneous medium 
	% e.g., tje conductivity and the Helmholtz equation. Some common characteristics:
	% 1. The domain is illuminated by a (point) source
	% 2. The solution can be represented using layer potentials
	% 3. Decomposition formula exists for the solution which allows to reconstruct the asymptotics (GPT).
	
	properties(SetAccess = protected)
		D = {}; % A list of inclusions (C2boundary)
		nbIncls = 0; % Number of inclusions
		
		cfg % acquisition configuration
	end
	
	methods
		function obj = Small_Inclusions(D, cfg)
			% D: a list of C2boundary objects
			% cfg: an acquisition configuration object
			
			if ~iscell(D)
				% if the input is a single inclusion
				obj = obj.addInclusion(D);
			else
				% or if the input is a list of inclusions
				for n=1:length(D)
					obj = obj.addInclusion(D{n});
				end
			end
			
			if ~isa(cfg, 'acq.mconfig')
				error('TypeError: must be an object of acq.mconfig');
			else
				obj.cfg = cfg;
			end
		end
		
		function obj = addInclusion(obj, D)
			% Add a small inclusion
			if ~isa(D, 'shape.C2boundary')
				error('Type error: the inclusion must be an object of C2boundary')
			end
			
			if obj.nbIncls >= 1
				if obj.D{obj.nbIncls}.nbPoints ~= D.nbPoints
					error('Value error: the number of boundary discretization points must be the same for all inclusions.');
				end
				
				if ~obj.check_inclusions(D)
					error('Inclusions must be separated from each other.');
				end
			end
			
			obj.D{obj.nbIncls+1} = D;
			obj.nbIncls = obj.nbIncls+1;
		end
		
		function val = check_inclusions(obj, D)
			val = 1;
			for n=1:obj.nbIncls
				val = val * D.isdisjoint(obj.D{n});
			end
		end
		
		function plot(obj, varargin)
			% plot first the inclusions
			for n=1:obj.nbIncls
				plot(obj.D{n}, varargin{:}); hold on;
			end
			% plot then the acquisition system
			plot(obj.cfg, varargin{:});
		end
	end
	
	methods(Abstract)
		% Simulate the perturbation of field u-U. u and U are fields with and without
		% inclusions. This perturbation is also the data measured by receptors. Remark that there
		% exists other type of perturbation, so forth other type of data.
		out = data_simulation(obj)
		
		% % Solve the forward problem by evaluating the solution using its
		% % representation. The meaning as well as the implementation of this function can differ
		% % according to the definition of the problem.
		% u = solve_forward()
	end
end
