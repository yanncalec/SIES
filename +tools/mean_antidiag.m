function antidiag=mean_antidiag(Z, ord)
% antidiag=mean_antidiag(Z, ord)
% Return the mean absolute value of the anti-diagonals 1..ord of
% matrix Z

antidiag = zeros(1,ord);
for at=1:ord
    npm = at+1;
    for n=1:at
        m = npm-n;
        antidiag(at) = antidiag(at) + abs(Z(n,m)); 
    end 
    antidiag(at) = antidiag(at)/at;
end 
