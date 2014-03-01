classdef mconfig     
% Abstract class for the geometric configuration of the acquisition system.
%  
% An acquisition system consists of sources and receivers. This class contains only the geometry
% information, eg, the position, the direction of the dipole source or the plane wave etc, but not
% the physical information, like the frequency. The reason is that the target to be measured is
% independent to the geometry of the acquisition system, but dependent to the frequency of the
% source, for example. The physical model of the source is included in the PDE.Small_Inclusions
% class.
%
% Sources and receivers are divided into groups. Inside a group all receivers are visible to all
% sources. No interaction exists between the receivers of group A and sources of group B. The
% number of sources (resp. receivers) are the same for all groups.
    
    properties(Access = protected)
        src_prv % coordinates of sources by group, a cell of array of dimension 2 X Ns
        rcv_prv % coordinates of receivers by group, a cell of array of dimension 2 X Nr
    end

    properties(SetAccess = protected)        
        %% Manually-set
        Ng % number of groups
        center % reference center of the mesurement system
        
        %% Auto-set
        Ns % number of sources in each group
        Nr % number of receivers  in each group
        Ns_total % total number of sources
        Nr_total % total number of receivers 
        data_dim % total dimension of the data
    end
    
    methods
        function val = get.Ns(obj)
            val = size(obj.src_prv{1},2);
        end
        
        function val = get.Nr(obj)
            val = size(obj.rcv_prv{1},2);
        end
        
        function val = get.Ns_total(obj)
            val = obj.Ns * obj.Ng;
        end
        
        function val = get.Nr_total(obj)
            val = obj.Nr * obj.Ng;
        end
        
        function val = get.data_dim(obj)
            val = obj.Ns*obj.Nr*obj.Ng;
        end

        %% Utility functions
        function [src, rcv] = group(obj, g)
        % Coordinates of the sources and receivers of the g-th group
            src = obj.src_prv{g};
            rcv = obj.rcv_prv{g};
        end
        
        function [gid, sid] = src_query(obj, s)
        % For the s-th source, get the group index and the index inside the group
            if s>obj.Ns_total || s<1
                error('Source index out of range');
            end
            gid = floor((s-1)/obj.Ns)+1;
            sid = s - (gid-1)*obj.Ns;
        end
        
        % function [gid, rid] = rcv_query(obj, r)
        % % For the r-th receiver, get the group index and the index inside the group
        %     if r>obj.Nr_total || r<1
        %         error('Receiver index out of range');
        %     end
        %     gid = floor((r-1)/obj.Nr)+1;
        %     rid = r - (gid-1)*obj.Nr;
        % end
        
        function val = src(obj,s)
        % Coordinates of the s-th source
            [gid, sid] = obj.src_query(s);
            
            toto = obj.src_prv{gid};
            val = toto(:,sid);
        end

        function val = rcv(obj, s)
        % Coordinates of receivers responding to the s-th source
            [gid, sid] = obj.src_query(s);
            
            val = obj.rcv_prv{gid};
        end

        function val = all_src(obj)
        % Coordinates of all sources, concatenated group by group
            val = [];
            for g=1:obj.Ng
                val = [val obj.src_prv{g}];
            end
        end

        function val = all_rcv(obj)
        % Coordinates of all receivers, concatenated group by group
            val = [];
            for g=1:obj.Ng
                val = [val obj.rcv_prv{g}];
            end
        end
        
        function plot(obj, varargin)
            for g=1:obj.Ng
                [src, rcv] = obj.group(g);
                plot(src(1,:), src(2,:), 'x'); hold on;
                plot(rcv(1,:), rcv(2,:), 'o');
            end
        end        

    end
    
end