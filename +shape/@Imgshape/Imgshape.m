classdef Imgshape < shape.C2boundary
% Class for importing shape from an image file
    
    methods
        function obj = Imgshape(fname, nbPoints, a, b)
        % Initialize a standard boundary from an image file
        % The non tied-off parameterization must be followed everywhere.
            if nargin < 4
                a=1; b=1;
            end
            
            im = shape.Imgshape.loadimgfnc(fname);                 % Load image
            
            [D0, theta0] = shape.Imgshape.boundarydet(0.5, im);    % Extract the boundary
            
            [points, tvec, avec, normal] = shape.C2boundary.rescale(D0, theta0, nbPoints, [a, b]);
            obj = obj@shape.C2boundary(points, tvec, avec, normal, []);

            % obj.points = points;
            % obj.nbPoints = nbPoints;
            % obj.tvec = tvec;
            % obj.avec = avec;
            % obj.normal = normal;
            % obj.center_of_mass = obj.get_center_of_mass();
            
            idx1 = strfind(fname,'/'); idxslash=idx1(end);
            idx2 = strfind(fname,'.'); idxdot=idx2(end);            
            obj.name_str = fname(idxslash+1:idxdot-1); % Generate the name string from filename                        
        end
    end    

    methods(Static)
        IM=loadimgfnc(FileName)
        [D,t] = boundarydet(acc, IM)
    end 
end
