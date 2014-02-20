classdef ROI
% Rectangular region of interest (ROI) for wavelet basis

    properties(SetAccess=protected)
        center % center of ROI
        width % size in x-axis
        height % size in y-axis

        suppx % support interval [a,b] in x-axis
        suppy % support interval [a,b] in y-axis
    end
    
    methods
        function val = get.suppx(obj)
            val=obj.center(1) + [-obj.width, obj.width]/2;
        end

        function val = get.suppy(obj)
            val=obj.center(2) + [-obj.height, obj.height]/2;
        end

        function obj = ROI(center, width, height)
            obj.center = center;
            obj.width = width;
            obj.height = height;            
        end 
    end
end