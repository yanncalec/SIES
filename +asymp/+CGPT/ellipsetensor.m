function exact = ellipsetensor(order,a,b,k)
% Calculate the exact CGPT of an ellipse.
% Inputs:
% -order: maximum order of CGPT
% -a,b: lengths of the semi-major and semi-minor axes of the ellipse
% -k: conductivity constant
% WARNING: the result with complex conductivity is not guaranteed, use
% theoretical_CGPT instead.

% ACTUAL ELLIPSE PT:
% We refer the reader to the PhD thesis of Mikyoung Lim.  In it she
% derives the exact calculation of the GPT for an ellipse. Those results
% are also reported in Ammari & Kang 2007 Book
%
%  Created 23-SEP-2011 by  Anton Bongio Karrman karrman@caltech.edu
%
exact = zeros(2*order,2*order) ;
for m = 1:order
	for n = m:order
		exact(2*m-1:2*m,2*n-1:2*n) = real(actualE(m,n,a,b,k)) ;
	end
end
%    exact = exact+(exact.*triu(ones(2*order,2*order),1))'     ;
exact = exact+(exact.*triu(ones(2*order,2*order),1))';
exact = exact * (k-1);
end

function matel = actualE(m,n,a,b,k)
p = atanh(b/a) ;
r = sqrt(a^2-b^2) ;
matel = zeros(2) ;
cm = constm(m,n,r) ;
if ((n-m)/2 >= 0) && (mod(m-n,2) == 0)
	caR = consta(m,n,p,k,0) ;
	caI = consta(m,n,p,k,1)  ;
end
if (mod(m,2) == 0) && (mod(n,2) == 0)
	matel(1,1) = cm*(evenSum(m,n,p,k,0)+caR) ;
	matel(2,2) = cm*(evenSum(m,n,p,k,1)+caI) ;
else if (mod(m,2) == 1) && (mod(n,2) == 1)
		matel(1,1) = cm*(oddSum(m,n,p,k,0)+caR) ;
		matel(2,2) = cm*(oddSum(m,n,p,k,1)+caI) ;
	end
end
end

function val = constm(m,n,r)
val = m*pi*r^(m+n)*2^(1-m-n) ;
end

function val = consta(m,n,p,k,til)
val = abs(sign(n-m+2))*((sign(n-m+2)+1)/2)...
	*nchoosek(n,(n-m)/2)*bf(m,p,k,til)*sinh(2*m*p) ;
end

function val = evenSum(m,n,p,k,til)
minval = min(n/2,(m-2)/2) ;
if (minval >= 1)
	t = 1:minval ;
	parta = 4*factorial(m-1)*factorial(n) ;
	parta = parta*bf(2*t,p,k,til).*t.*sinh(4*t*p) ;
	partb = (m+2*t).*factorial(m/2-t).*factorial(m/2-1+t) ;
	partb = partb.*factorial(n/2+t).*factorial(n/2-t) ;
	val = sum(parta./partb) ;
else
	val = 0 ;
end
end

function val = oddSum(m,n,p,k,til)
minval = min((n-1)/2,(m-3)/2) ;
if (minval >= 0)
	t = 0:minval ;
	parta = factorial(m-1)*factorial(n) ;
	parta = parta*bf(2*t+1,p,k,til).*(4*t+2).*sinh((4*t+2)*p) ;
	partb = (m+2*t+1).*factorial((m-1)/2-t).*factorial((m-1)/2+t) ;
	partb = partb.*factorial((n+1)/2+t).*factorial((n-1)/2-t) ;
	val = sum(parta./partb) ;
else
	val = 0 ;
end
end

function val = bf(v,p,k,til)
if (til == 0)
	val = (sinh(v*p)+cosh(v*p))./(k*sinh(v*p)+cosh(v*p)) ;
else
	val = (sinh(v*p)+cosh(v*p))./(sinh(v*p)+k*cosh(v*p)) ;
end
end
%
%  Created 23-SEP-2011 by  Anton Bongio Karrman karrman@caltech.edu
%

