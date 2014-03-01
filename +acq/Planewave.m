classdef Planewave < acq.mconfig
    % Configuration with fixed receivers and plane wave sources.
    %
    % All receivers (resp. sources) are equally distributed on a circle (resp. between [0,2pi]).  All
    % receivers are visible to every source and are independent to sources. The sources are identified
    % by their directions.
        
    properties(SetAccess = protected)
        radius_src = 1% radius of measurement circle for sources
        
        radius_rcv % radius of measurement circle for receivers
        equispaced = 0 % 1 if sources and receivers are equispaced
    end


    methods        
        function obj = Planewave(Z, Rr, Nr, Ns, viewmode, grouped)
            % INPUTS:
            % Ns_total: total number of sources
            % Z: center of the measurement circle
            % Rr: radius of the circle receivers
            % Nr_total: total number of receivers
            % viewmode: three dimension vectors [Na, theta, aov] specifying
            %        Na: sources are divided into Na number of separated arcs
            %        theta: angular aperture of each arc
            %        aov: angle covered by the starting point of all arcs
            
            if nargin < 6
                grouped = 0;
            end
            
            if nargin < 5
                viewmode = [1, 2*pi, 2*pi];
            end
            
            Na = viewmode(1); theta = viewmode(2); aov = viewmode(3);
            
            [Xs, ~, Xscell] = acq.src_rcv_circle(Na, Ns, 1, [0,0]', theta, aov);
            [Xr, ~, Xrcell] = acq.src_rcv_circle(Na, Nr, Rr, Z, theta, aov);
            
            obj.center = Z;
            obj.radius_rcv = Rr;
            
            if grouped
                obj.Ng = Na;
                obj.src_prv = Xscell;
                obj.rcv_prv = Xrcell;
            else            
                obj.Ng = 1;            
                obj.src_prv = {Xs};
                obj.rcv_prv = {Xr};
            end                
            
            if obj.Ng==1 && theta==2*pi
                obj.equispaced = 1;
            end

        end
        
        function plot(obj, varargin)
            % The plot for source is marked by 'x' and indicates only their direction between
            % [0,2*pi].
            
            for g=1:obj.Ng
                [src, rcv] = obj.group(g);

                src = repmat(obj.center, 1, obj.Ns) + obj.radius_rcv * src;
                plot(src(1,:), src(2,:), 'x'); hold on;
                
                rcv = obj.radius_src * rcv;
                plot(rcv(1,:), rcv(2,:), 'o');
            end
            plot(obj.center(1), obj.center(2), 'r*');
        end
    end
end

