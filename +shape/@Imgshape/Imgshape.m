classdef Imgshape < shape.C2boundary
% Class for importing shape from an image file
    
    methods
        function obj = Imgshape(fname, nbPoints, a, b, hwidth)
        % Initialize a standard boundary from an image file
        % The non tied-off parameterization must be followed everywhere.

            if nargin<5
                hwidth = 10;
            end

            if nargin < 4
                a=1; b=1;
            end
            
            im = shape.Imgshape.loadimgfnc(fname);                 % Load image
            
            [D0, theta0] = shape.Imgshape.boundarydet(0.5, im);    % Extract the boundary
            
            [points, ~] = shape.C2boundary.rescale(D0, theta0, nbPoints, [a, b]);
            
            [points, tvec, avec, normal] = shape.C2boundary.smooth_out_singularity(points, [0,0]', hwidth, [a,b]);
            
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