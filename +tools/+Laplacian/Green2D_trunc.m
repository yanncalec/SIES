function [G] = Green2D_trunc(X, Y, epsilon)
    
    if nargin < 3
        epsilon = 1e-7;
    end
    
    if size(X,1)~=2 || size(Y,1)~=2
        error('The inputs X and Y must have two rows!');
    end
    % X = reshape(X, 2, []); 
    % Y = reshape(Y, 2, []);

    XY1 = tools.tensorplus(X(1,:), -Y(1,:));
    XY2 = tools.tensorplus(X(2,:), -Y(2,:));

    DN = sqrt(XY1.^2 + XY2.^2);
    idx = find(DN<epsilon);
    G = 1/2/pi*log(DN);
    G(idx) = 1/2/pi*log(epsilon);
end
