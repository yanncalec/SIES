classdef Flower < shape.C2boundary
    % Class for flowers (rotational symmetric shapes) based on an ellipse
    
    properties
        nbPetals % number of petals
        axis_a % length of the semi-major axis
        axis_b % length of the semi-major axis
        phi % orientation of the ellipse
        epsilon % size of pertubation
        pertb % the exponent used in the perturbation

        % symm % number of rotational symmetry
        tau % percentage of damage
    end
    
    methods
        function obj = Flower(a, b, nbPoints, nbPetals, epsilon, tau)
            if nargin<6
                tau = 0;
            end
            if nargin < 5
                epsilon = 0.3;
            end
            if nargin<4
                nbPetals = 5;
            end            
            
            x0 = 0; y0 = 0; phi = 0;
            
            com = [0,0]';
            pertb = 1;
            
            if tau == 0
                [~,points,tvec,~,normal,avec,sigma,~] = shape.Flower.make_flower(x0,y0,nbPetals, pertb, a, b, phi, epsilon, nbPoints);
            elseif tau > 0 && tau < 1
                [~,points,tvec,~,normal,avec,sigma] = shape.Flower.make_damaged_flower(x0,y0,nbPetals,a,b,phi,epsilon,nbPoints,tau) ;

                com = shape.C2boundary.get_com(points, tvec, normal);
            else
                error('Value error: the percentage of damage must be between 0 and 1.')
            end            
            
            obj = obj@shape.C2boundary(points, tvec, avec, normal, com, 'Flower');

            obj.nbPetals = nbPetals;
            obj.axis_a = a; obj.axis_b = b;
            obj.phi = phi; 
            obj.epsilon = epsilon; obj.tau = tau;
        end
        
        function obj = mtimes(obj, s)
        % Overload of the operator *
            obj = mtimes@shape.C2boundary(obj, s);
            obj.axis_a = obj.axis_a * s;
            obj.axis_b = obj.axis_b * s;            
        end                
        
        function obj = lt(obj, phi)
        % Redefine the < operator as the rotation of the boundary
            obj = lt@shape.C2boundary(obj,phi);
            obj.phi = obj.phi+phi;
        end
        
    end
    
    methods(Static)
        [theta,D,tvec,norm_Tan_square,normal,avec,Sigma,Ksymm] = make_flower(x0,y0,n,k,a,b,phi, ...
                                                          epsilon,N)
        
        [theta,D,tvec,norm_Tan_square,normal,avec,Sigma] = make_damaged_flower(x0,y0,n,a,b,phi, ...
                                                          epsilon,N,tau)        
    end
end

