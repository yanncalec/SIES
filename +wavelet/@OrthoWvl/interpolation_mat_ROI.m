function [P, dPX, dPY] = interpolation_mat_ROI(obj, X, nbScl)
% Make wavelet interpolation matrix on a ROI
% INPUTS:
% X: coordinates of interpolation, 2-by-N
% J: Jmax..Jmin-1, scale of interpolation. The scaling functions of this scale
% are used as interpolation functions.
% OUTPUTS:
% Phi: interpolation matrix, N-by-M, with M the number of active scaling functions
% dPhiX, dPhiY:  interpolation matrix of derivatives

    [Phi0, dPhiX0, dPhiY0] = wavelet.OrthoWvl.evaluation(X, obj.Jmax, ...
                                                      obj.ApprxSpace{1}.pmeshpoints, ...
                                                      obj.Tphi, obj.Tphi, obj.Gphi, ...
                                                      obj.Gphi);
    [M,N]=size(Phi0);
    Phi=cell(nbScl, 3);
    dPhiX=cell(nbScl, 3);
    dPhiY=cell(nbScl, 3);
    
    for j=1:nbScl
        J = obj.Jmax-j+1;
        [Phi{j,1}, dPhiX{j,1}, dPhiY{j,1}] = wavelet.OrthoWvl.evaluation(X, ...
                                                          J,obj.DetailSpace{j, ...
                            1}.pmeshpoints, obj.Tphi, obj.Tpsi, obj.Gphi, ...
                                                          obj.Gpsi);
        [Phi{j,2}, dPhiX{j,2}, dPhiY{j,2}] = wavelet.OrthoWvl.evaluation(X, ...
                                                          J,obj.DetailSpace{j, ...
                            2}.pmeshpoints, obj.Tphi, obj.Tpsi, obj.Gphi, ...
                                                          obj.Gpsi);
        [Phi{j,3}, dPhiX{j,3}, dPhiY{j,3}] = wavelet.OrthoWvl.evaluation(X, ...
                                                          J,obj.DetailSpace{j, ...
                            3}.pmeshpoints, obj.Tphi, obj.Tpsi, obj.Gphi, ...
                                                          obj.Gpsi);

        N = N+size(Phi{j,1},2)+size(Phi{j,2},2)+size(Phi{j,3},2);
    end                                                              

    P=zeros(M,N);
    dPX=zeros(M,N);
    dPY=zeros(M,N);

    P(:,1:size(Phi0,2)) = Phi0;    
    dPX(:,1:size(Phi0,2)) = dPhiX0;
    dPY(:,1:size(Phi0,2)) = dPhiY0;
    
    ptr = size(Phi0,2)+1;
    
    for j=1:nbScl
        for k=1:3
            P(:, ptr:ptr+size(Phi{j,k},2)-1) = Phi{j,k};
            dPY(:, ptr:ptr+size(Phi{j,k},2)-1) = dPhiX{j,k};
            dPY(:, ptr:ptr+size(Phi{j,k},2)-1) = dPhiY{j,k};
            ptr = ptr+size(Phi{j,k},2);
        end
    end
end
