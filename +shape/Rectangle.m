classdef Rectangle < shape.C2boundary
    % Class for rectangular shape
    properties(SetAccess = protected)
        width
        height
    end
    
    methods
        function obj = Rectangle(a,b,nbPoints)
        % This function creates a structure representing an anomaly which has the
        % shape of a rectangle.
        % INPUT : - a = height,b = width
        %         - nPoints = number of discretization points
            
            t1 = b/(a+b)/2; t2 = a/(a+b)/2; t3 = t1; t4 = t2;
            n1 = floor(t1*nbPoints); n2 = floor(t2*nbPoints); 
            n3 = floor(t3*nbPoints); n4 = nbPoints-n1-n2-n3;
            
            A = [-b ; a]/2; 
            B = [-b; -a]/2;
            C = [b ; -a]/2; 
            D = [b ; a]/2; 

            t = (0:n1-1)/n1;
            AB = repmat(A,1,n1) + repmat(B-A,1,n1) .* repmat(t,2,1);
            t = (0:n2-1)/n2;
            BC = repmat(B,1,n2) + repmat(C-B,1,n2) .* repmat(t,2,1);
            t = (0:n3-1)/n3;
            CD = repmat(C,1,n3) + repmat(D-C,1,n3) .* repmat(t,2,1);
            t = (0:n4-1)/n4;
            DA = repmat(D,1,n4) + repmat(A-D,1,n4) .* repmat(t,2,1);

            points = [AB BC CD DA] ;
            tvec = [repmat((B-A),1,n1)/t1 repmat((C-B),1,n2)/t2 repmat((D-C),1,n3)/t3 repmat((A-D),1,n4)/t4]/2/pi ; % velocity vector

            rotation = [[0 1];[-1 0]] ;
            normal = rotation*tvec ;
            normal = normal./repmat(sqrt(normal(1,:).^2+normal(2,:).^2),2,1) ; % normal vector

            avec = zeros(2,nbPoints) ; % acceleration vector
            com = [0, 0]';

            if a==b
                name_str = 'Square';
            else
                name_str = 'Rectangle';
            end
            obj = obj@shape.C2boundary(points, tvec, avec, normal, com, name_str);
            
            obj.width = b;
            obj.height = a;
        end        
        
        function obj = mtimes(obj, s)
        % Overload of the operator *
            obj = mtimes@shape.C2boundary(obj, s);
            obj.width = obj.width * s;
            obj.height = obj.height * s;            
        end                
        
    end    
end

    
% Old version by Thomas:
%     function D = Rectangle(a,b,center,phi,M)
%         % This function creates a structure representing an anomaly which has the
%         % shape of a rectangle.
%         % INPUT : - a,b = length of the border
%     		%         - center = coordinates of the center
%         %         - phi = rotation angle
%         %         - M = number of discretization points

%         x0 = center(1); y0 = center(2);            

%         rot = [[cos(phi) -sin(phi)];[sin(phi) cos(phi)]] ;

%         D.nbPoints = 2*M ;
%         Nx = floor(a/(a+b)*M) ;
%         Ny = M-Nx ;

%         dSigmaX = a/Nx ;
%         westX = (-a/2)*ones(1,Ny) ;
%         northX = linspace(-a/2+dSigmaX,+a/2-dSigmaX,Nx) ;
%         eastX = (+a/2)*ones(1,Ny) ;
%         southX = linspace(+a/2-dSigmaX,-a/2+dSigmaX,Nx) ;

%         dSigmaY = b/Ny ;
%         westY = linspace(-b/2+dSigmaY,+b/2-dSigmaY,Ny) ;
%         northY = (+b/2)*ones(1,Nx) ;
%         eastY = linspace(+b/2-dSigmaY,-b/2+dSigmaY,Ny) ;
%         southY = (-b/2)*ones(1,Nx) ;

%         D.points = repmat([x0 ; y0],1,2*M) + ...
%             rot*[westX northX eastX southX ; ...
%             westY northY eastY southY ] ;

%         D.tvec = rot*[zeros(1,Ny) dSigmaX*ones(1,Nx) zeros(1,Ny) -dSigmaX*ones(1,Nx);...
%             dSigmaY*ones(1,Ny) zeros(1,Nx) -dSigmaY*ones(1,Ny) zeros(1,Nx)] ; % velocity vector
%         %D.norm_Tan_square = D.tvec(1,:).^2 + D.tvec(2,:).^2 ;
%         rotation = [[0 -1];[1 0]] ;
%         normal = rotation*D.tvec ;
%         D.normal = normal./repmat(sqrt(normal(1,:).^2+normal(2,:).^2),2,1) ;
%         D.avec = zeros(2,2*M) ; % acceleration vector
%         %D.sigma = [dSigmaY*ones(1,Ny) dSigmaX*ones(1,Nx) dSigmaY*ones(1,Ny) dSigmaX*ones(1,Nx)] ;

%         D.center_of_mass = center(:);

%         if a==b
%             D.name_str = 'Square';
%         else
%             D.name_str = 'Rectangle';
%         end
%     end        
% end
