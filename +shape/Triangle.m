classdef Triangle < shape.C2boundary
    % Isosceles triangle (with smoothed corner)
    
    properties
        lside % length of the equal side
        angl % angle between the two equal sides
    end
    
    methods
        function obj = Triangle(a, angl, nbPoints)
        % INPUT :
        % a: length of the equal side
        % angl: angle between the two equal sides
        % nbPoints: number of discretization points

            h = a*cos(angl/2); % height of the triangle
            b = a*sin(angl/2); % 2b is the length of the 3rd side
            
            t1 = a/(a+b)/2; t2 = b/(a+b); t3 = t1;
            n1 = floor(t1*nbPoints); n2 = floor(t2*nbPoints); n3 = nbPoints-n1-n2;
            
            A = [ 0 ; 2/3*h ] ; B = [ - b ; -h/3 ] ; C = [ b ; -h/3 ] ;

            t = (0:n1-1)/n1;
            AB = repmat(A,1,n1) + repmat(B-A,1,n1) .* repmat(t,2,1);
            t = (0:n2-1)/n2;
            BC = repmat(B,1,n2) + repmat(C-B,1,n2) .* repmat(t,2,1);
            t = (0:n3-1)/n3;
            CA = repmat(C,1,n3) + repmat(A-C,1,n3) .* repmat(t,2,1);
            % CA = (1-t)*repmat(C,1,M) + t*repmat(A,1,M);

            points = [AB BC CA] ;
            tvec = [repmat((B-A),1,n1)/t1 repmat((C-B),1,n2)/t2 repmat((A-C),1,n3)/t3]/2/pi ; % velocity vector

            rotation = [[0 1];[-1 0]] ;
            normal = rotation*tvec ;

            normal = normal./repmat(sqrt(normal(1,:).^2+normal(2,:).^2),2,1) ;
            avec = zeros(2,nbPoints) ; % acceleration vector
            com = [0, 0]'; % the triangle is centered at the origine
            
            obj = obj@shape.C2boundary(points, tvec, avec, normal, com, 'Triangle');
            obj.lside = a;
            obj.angl = angl;
        end
    end
end

