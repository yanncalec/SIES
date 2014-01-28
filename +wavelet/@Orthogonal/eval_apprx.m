function val = eval_apprx(obj, X)
% Evaluate scaling functions of the approximation space. The scale number for the
% approximation space is Jmax. Only those wavelets who contribute to the ROI
% are active and their position number n define the function phi_{Jmax, n}
% which is the tensor-product of two 1D scaling functions.
    
    N =  obj.ApprxSpace{1}.pmeshpoints;
    val = wavelet.OrthoWvl.evaluation(X, obj.Jmax, N, obj.Tphi, obj.Tphi, obj.Gphi, obj.Gphi);
end
