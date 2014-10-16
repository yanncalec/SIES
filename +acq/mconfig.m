classdef mconfig
	% Abstract class for the configuration of acquisition system.
	%
	% An acquisition system consists of sources and receivers. This class
	% contains only the geometrical properties, eg, the position, the direction of the dipole source 
	% or of the plane wave etc, but contains no physical properties like the frequency of the source. 
	% The reason of this design is that the target to be measured is independent to the geometry of 
	% the acquisition system, but may dependent to the frequency. The physical
	% properties are included in the PDE.Small_Inclusions class.
	%
	% Sources and receivers are divided into groups. Inside a group all sources and receivers are
	% mutually visible. No interaction exists between the receivers of one group and sources of 
	% another group. The number of sources (resp. receivers) are the same for all groups.
	
	properties(Access = protected)
		src_prv % coordinates of sources by group, a cell of array of dimension 2 X Ns
		rcv_prv % coordinates of receivers by group, a cell of array of dimension 2 X Nr
	end
	
	properties(SetAccess = protected)
		%% Manually set
		Ng % number of groups
		center % reference center of the mesurement system
		
		%% Automatically set
		Ns % number of sources in each group
		Nr % number of receivers in each group
		Ns_total % total number of sources
		Nr_total % total number of receivers
		data_dim % total dimension of the data matrix
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
			val = obj.Ns*obj.Ng*obj.Nr; % There are in total Ns*Ng sources, and to each source are associated Nr receivers (those in the same group).
		end
		
		%% Utility functions
		function [src, rcv] = group(obj, g)
			% Coordinates of the sources and receivers of the g-th group
			src = obj.src_prv{g};
			rcv = obj.rcv_prv{g};
		end
		
		function [gid, sid] = src_query(obj, s)
			% For the s-th source, get the group index and the index inside the group
			if length(s)~=1 || s>obj.Ns_total || s<1
				error('Non-scalar source index or source index out of range');
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
		
		function val = src(obj,sidx)
			% Coordinates of the s-th source
			val = zeros(2,length(sidx));
			
			for s=1:length(sidx)
				[gid, sid] = obj.src_query(sidx(s));
				toto = obj.src_prv{gid};
				val(:,s) = toto(:,sid);
			end
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