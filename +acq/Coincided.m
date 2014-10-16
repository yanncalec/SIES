classdef Coincided < acq.Concentric
	% Coincided sources/receivers distributed on a circle.
	%
	% Each source acts also as a receiver. The sources are equally spaced inside the arcs which are 
	% uniformly distributed on a circle. This is the simplest configuration for the conductivity problem.
	
	methods
		function obj = Coincided(Z, Rs, Ns, viewmode, grouped, neutCoeff, neutRad)
			% INPUTS:
			% -Z: center of the measurement circle
			% -Rs: radius of the measurement circle
			% -Ns: number of sources/receivers per arc
			% -viewmode...neutRad: same as in the class Concentric.m
			
			if nargin < 7
				neutRad = 0.01;
			end
			
			if nargin < 6
				neutCoeff = [];
			end
			
			if nargin < 5
				grouped = 0;
			end
			
			if nargin < 4
				viewmode = [1, 2*pi, 2*pi];
			end
			
			obj = obj@acq.Concentric(Z, Rs, Ns, Rs, Ns, viewmode, grouped, neutCoeff, neutRad);
		end
	end
end
