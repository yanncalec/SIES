function Y = mdiff(X,d)
% Y = mdiff(X,d)
% Middle point difference of a matrix X in direction d(d=1 for row or 2 for column).

if nargin<2
    d=1;
end

Y=zeros(size(X));
[M,N] = size(X);

if M>1 && N>1    
    if d==1
        Y(2:end-1,:) = (X(3:end,:)-X(1:end-2,:))/2;
        Y(1,:) = X(2,:)/2;
        Y(end,:) = -X(end-1,:)/2;
    else
        Y(:,2:end-1,:) = (X(:,3:end)-X(:,1:end-2))/2;
        Y(:,1) = X(:,2)/2;
        Y(:,end) = -X(:,end-1)/2;        
    end
else
    Y(2:end-1) = (X(3:end)-X(1:end-2))/2;
    Y(1) = X(2)/2;
    Y(end) = -X(end-1)/2;
end