%% Related to the hankel function

%% H
ord=0:20;
z=0.1+0.1*1i;
hankelf = besselj(ord, z) + 1i*bessely(ord, z);

figure; semilogy(abs(hankelf));


%% Verify the inequality of 
% sum_(m>k) (c/m)^m <= (c/k)^k/(1+ln(k/c))

for nn=1:100
    c = rand*10;
    k = ceil(c*exp(-1));
    % k = max(c*exp(-1), ceil(c*exp(-1)));
    N=2000;

    k0=k+1;
    M=k0:(k0+N);
    ss = (c./M).^M;
    s=sum(ss);
    res(nn)=(c/k)^k/(1+log(k/c)) - s;
    %res(nn)=(c/k)^k - s;
    % res(nn)=exp(k) - s;
end

min(res)

