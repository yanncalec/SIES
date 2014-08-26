classdef Imgshape < shape.C2boundary
    % Class for importing shape from an image file
    
    methods
        function obj = Imgshape(fname, nbPoints, a, b, dspl)
            % Initialize a standard boundary from an image file
            % Inputs:
            % fname: name of the image file
            % nbPoints: number of boundary points
            % a,b: width and height of the final shape
            % dspl: down-sampling factor for smoothing the corners
            
            if nargin < 5
                dspl = 5;
            end
            
            if nargin < 4
                a=1; b=1;
            end
            
            im = shape.Imgshape.loadimgfnc(fname);                 % Load image
            
            [points0, theta0] = shape.Imgshape.boundarydet(0.5, im);    % Extract the boundary
            
            if dspl >= 1
                [points, tvec, avec, normal] = shape.C2boundary.rescale(points0, theta0, nbPoints, [a, b], dspl);
            else
                [points, tvec, avec, normal] = shape.C2boundary.rescale_diff(points0, theta0, nbPoints, [a, b]);
            end
        
            idx1 = strfind(fname,'/'); idxslash=idx1(end);
            idx2 = strfind(fname,'.'); idxdot=idx2(end);
            name_str = fname(idxslash+1:idxdot-1); % Generate the name string from filename
            
            obj = obj@shape.C2boundary(points, tvec, avec, normal, [], name_str);
        end
    end
    
    methods(Static)
        IM=loadimgfnc(FileName)
        [D,t] = boundarydet(acc, IM)
    end
end