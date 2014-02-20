function val = eval_detail(obj, X, J, k)
% Evaluate wavelet functions of the detail space of scale number J. Only those
% wavelets who contribute to the ROI are active and their position number n
% define the function psi_{J, n} which is the tensor-product of two 1D
% wavelet functions.

    if J<obj.Jmin || J>obj.Jmax
        error('Scale number out of valid range!');
    end

    j = obj.Jmax - J + 1;
    N = obj.DetailSpace{j,k}.pmeshpoints;
    % scl = obj.Jmax-j+1;

    if k==1
        val = wavelet.OrthoWvl.evaluation(X, J, N, obj.Tphi, obj.Tpsi, ...
                                          obj.Gphi, obj.Gpsi);
    elseif k==2
        val = wavelet.OrthoWvl.evaluation(X, J, N, obj.Tpsi, obj.Tphi, ...
                                          obj.Gpsi, obj.Gphi);
    else                
        val = wavelet.OrthoWvl.evaluation(X, J, N, obj.Tpsi, obj.Tpsi, ...
                                          obj.Gpsi, obj.Gpsi);
    end                
end
