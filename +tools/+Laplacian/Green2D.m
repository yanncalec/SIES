function [G] = Green2D(X, Y)
% [G] = Green2D(X, Y)
%
% 2D Green function G(x) = 1/2/pi * log|x| evaluated at the points X-Y. The inputs X,
% Y are arrays of dimension 2X?, and the output G is a matrix whose (m,n) term is G
% evlauated at (X(:,m) - Y(:,n)).

    if size(X,1)~=2 || size(Y,1)~=2
        error('The inputs X and Y must have two rows!');
    end
    % X = reshape(X, 2, []); 
    % Y = reshape(Y, 2, []);

    XY1 = tools.tensorplus(X(1,:), -Y(1,:));
    XY2 = tools.tensorplus(X(2,:), -Y(2,:));

    DN = XY1.^2 + XY2.^2;
    G = 1/4/pi*log(DN);

end
