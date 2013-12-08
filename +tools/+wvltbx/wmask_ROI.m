function [M, Mc] = wmask_ROI(wname, N0, N1, nlvl)

    [~, L0] = wavedec(ones(1, N0), nlvl, wname);
    [~, L1] = wavedec(ones(1, N1), nlvl, wname);

    P=cumsum(L0(1:end-1))+1;
    cidx = wkeep(1:L0(1), L1(1));
    mask = zeros(1, L0(1)); mask(cidx)=1;

    M=[]; M{1}=mask;

    for n = 1:length(P)-1
        cidx = wkeep(1:L0(n+1), L1(n+1));
        mask = zeros(1, L0(n+1)); mask(cidx)=1;
        
        M{n+1}=mask;
    end

    Mc = M;
    M=cell2mat(Mc);
end