function [hA,hB]=hp_coeff(ordmin, ordmax)

hA = zeros(ordmax-ordmin+1);
hB = zeros(ordmax-ordmin+1);

for a1=ordmin:ordmax
    row = a1-ordmin+1;
    for a2=ordmin:ordmax
        col = a2-ordmin+1;
        
        [A, B] = GPT.get_hhp_coeff(a1+a2);
        hA(row, col) = A(col+1);
        hB(row, col) = B(col+1);    
    end 
end
