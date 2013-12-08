function M = disktensor(order,r,k)
% [exact,val1] = disktensor(order,r,k)
% INPUTS:
% order: maximum order of the CGPT
% r: radius of the disk
% k: conductivity constant
% OUTPUT:
% M: CGPT matrix
%
% We refer the reader to the PhD thesis of Mikyoung Lim.  In it she
% derives the exact calculation of the GPT for a disk:
%  
%  Created 23-SEP-2011 by  Anton Bongio Karrman karrman@caltech.edu
%

    n = (1:order) ;
    %    val1 = pi*n.*r.^(2*n)*2/(k+1)
    val1 = pi*n.*r.^(2*n)*2*(k-1)/(k+1) ;
    val2 = [val1; val1] ;
    M = diag(val2(:)) ;
end

%  
%  Created 23-SEP-2011 by  Anton Bongio Karrman karrman@caltech.edu
%

