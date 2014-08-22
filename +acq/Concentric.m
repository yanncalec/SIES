classdef Concentric < acq.mconfig
    % Abstract class for concentric configuration. All sources (resp. receivers) are placed on a
    % circle. The circle of sources and receivers are concentric. Sources are eventually divided into
    % groups, so that there is no communication between different groups, ie, the receivers of one group
    % cannot listen to sources of another group.
    
    properties(SetAccess = protected)
        radius_src % radius of measurement circle for sources
        radius_rcv % radius of measurement circle for receivers
        equispaced = 0 % 1 if sources and receivers are equispaced
        
        neutCoeff % Coefficient {a_j}_j of the neutrality condition
        nbDirac % Neutrality condition f(x) = sum_j=1^nbDirac a_j dirac(x-x_{s,j}), with sum_j a_j = 0
        neutRad = 0.1 % the positions of dirac x_{s,j} are distributed on a tangent segment at x_s of length proportional to neutRad
    end
    
    methods
        function obj = Concentric(Z, Rs, Ns, Rr, Nr, viewmode, grouped, neutCoeff, neutRad)
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
            % neutCoeff: coefficients of the neutrality condition
            % neutRad: a smalle positive number which determines the
            % distance between Diracs of the neutrality source:
            % neutRad*radius_src. Default value: 0.01
            
            if nargin < 9
                obj.neutRad = 0.01;
            end
            
            if nargin < 8 || length(neutCoeff) <= 1
                obj.neutCoeff = 1; % No neutrality condition in case of one Dirac source
            else
                if sum(neutCoeff) ~= 0 || max(abs(neutCoeff)) == 0
                    error('Coefficients of Diracs must be non zero and satisfy the neutrality condition (sum=0)!')
                end
                obj.neutCoeff = neutCoeff;
            end
            
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
        
        function val = get.nbDirac(obj)
            val = length(obj.neutCoeff);
        end
        
        function val = neutSrc(obj, s)
            % Get the positions (Diracs) of the s-th source fulfilling the neutrality
            % condition. The Diracs are distributed on a segment of length
            % obj.neutRad * obj.radius_src in the tangent direction to the
            % source circle.
            
            psrc = obj.src(s);
            
            if obj.nbDirac == 1
                val = psrc;
            else
                val = zeros(2, obj.nbDirac);
                L = obj.radius_src * obj.neutRad; % length of the segment
                toto = psrc - obj.center;
                q = [toto(2), -toto(1)]'; q = q/norm(q) * L;
                
                for n = 1:obj.nbDirac
                    val(:,n) = psrc + (n-1)/obj.nbDirac * q;
                end
            end
            %             tt = exp(1i* (0:obj.nbDirac-1)/obj.nbDirac * 2*pi);
            %             zz = obj.neutRad * obj.radius_src * tt;
            %             val = [real(zz(didx)); imag(zz(didx))] + psrc;
        end        
    end
end
