classdef Concentric < acq.mconfig
% Abstract class for concentric configuration. All sources (resp. receivers) are placed on a
% circle. The circle of sources and receivers are concentric. Sources are eventually divided into
% groups, so that there is no communication between different groups, ie, the receivers of one group
% cannot listen to sources of another group.

    properties(SetAccess = protected)
        radius_src % radius of measurement circle for sources
        radius_rcv % radius of measurement circle for receivers
        equispaced = 0 % 1 if sources and receivers are equispaced
    end

    methods
        function obj = Concentric(Z, Rs, Ns, Rr, Nr, viewmode, grouped)
        % INPUTS:
        % Z: center of the measurement circle
        % Rs, Rr: radius of the source/receiver measurement circle
        % Ns, Nr: number of source/receiver per arc
        % viewmode: three dimension vectors [Na, theta, aov] for sources/receivers specifying
        %        Na: s/r are divided into Na number of separated arcs
        %        theta: angular aperture of each arc
        %        aov: angle covered by the starting point of all arcs
        % grouped: if true each arc will be treated as an independant group (limited view), if false the whole
        % set of receivers will be visible for all sources (full view but sparse array).
            
            if nargin < 7
                grouped = 0;
            end

            if nargin < 6
                viewmode = [1, 2*pi, 2*pi];
            end

            Na = viewmode(1); theta = viewmode(2); aov = viewmode(3);
            
            [Xs, ~, Xscell] = acq.src_rcv_circle(Na, Ns, Rs, Z, theta, aov);
            [Xr, ~, Xrcell] = acq.src_rcv_circle(Na, Nr, Rr, Z, theta, aov);
            
            obj.center = Z;
            obj.radius_src = Rs;
            obj.radius_rcv = Rr;
            
            if grouped
                % obj.Ns = Ns;
                % obj.Nr = Nr;
                obj.Ng = Na;
                obj.src_prv = Xscell;
                obj.rcv_prv = Xrcell;
            else            
                % obj.Ns = Ns*Na;
                % obj.Nr = Nr*Na;
                obj.Ng = 1;            
                obj.src_prv = {Xs};
                obj.rcv_prv = {Xr};
            end                
            
            if obj.Ng==1 && theta==2*pi
                obj.equispaced = 1;
            end
        end
        
        function plot(obj, varargin)            
            plot@acq.mconfig(obj, varargin);
            plot(obj.center(1), obj.center(2), 'r*');
        end        
    end
end
