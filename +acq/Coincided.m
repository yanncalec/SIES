classdef Coincided < acq.Concentric
% Coincided sources/receivers distributed on a circle.
%
% Each source acts also as receivers inside each group. The sources are equally spaced as pieces of
% arcs on a circle of measurement. This is the simplest configuration for the conductivity problem.
    
    methods
        function obj = Coincided(Z, Rs, Ns, viewmode, grouped)
        % INPUTS:
        % Z: center of the measurement circle
        % Rs: radius of the measurement circle
        % Ns: number of sources/receivers per arc
        % viewmode: three dimension vectors [Na, theta, aov] specifying
        %        Na: s/r are divided into Na number of separated arcs
        %        theta: angular aperture of each arc
        %        aov: angle covered by the starting point of all arcs
        % grouped: if true each arc will be treated as an independant group (limited view), if false the whole
        % set of receivers will be visible for all sources (full view but sparse array).
            
            if nargin < 5
                grouped = 0;
            end

            if nargin < 4
                viewmode = [1, 2*pi, 2*pi];
            end

            obj = obj@acq.Concentric(Z, Rs, Ns, Rs, Ns, viewmode, grouped);
        end        
    end
end
