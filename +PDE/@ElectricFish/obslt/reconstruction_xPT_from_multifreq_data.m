function out = reconstruction_xPT_from_multifreq_data(MSR, Cur, PP_SFR, Cur_bg, cfg, ord, maxiter, tol, symmode)
% out = reconstruction_xPT_from_multifreq_data(MSR, Cur, PP_SFR, fpsi_bg, cfg, ord, maxiter, tol, symmode)
% Reconstruct CGPTs, PT and imag(PT) from data. 
%
% Inputs:
% MSR, Cur, PP_SFR: noisy data, output of the function data_simulation of an ElectricFish object
% fpsi_bg: background field, output of the function data_simulation of an ElectricFish object
% cfg: acquisition configuration
% ord: order of CGPT reconstruction
% maxiter: maximum number of iteration for lsqr iteration
% tol: tolerance of lsqr algorithm
% symmode: force the solution to be symmetric

%% Reconstruction of CGPT and PT for all frequencies

for n = 1:length(MSR)
    % %% Reconstruction of high order (>=2) CGPT and identification by matching the shape-descriptors
    
    % Reconstructed value
    % Don't forget to take the transpose of the data, each row needs to correspond to a source
    out1 = PDE.ElectricFish.reconstruct_CGPT(cfg, MSR{n}, Cur{n}, ord, maxiter, tol, symmode);
    out.CGPT{n} = out1.CGPT;
    
    % %% Reconstruction of PT from the post-processed SFR (PP_SFR) matrix and identification by matching the Polarization Tensors
    % This is more robust to noise than the CGPT reconstruction?
    
    out2 = PDE.ElectricFish.reconstruct_PT(cfg, PP_SFR{n}, fpsi_bg, maxiter, tol, symmode);
    out.PT{n} = out2.PT;
    
    % %% Reconstruction of the imaginary part of PT
    out3 = PDE.ElectricFish.reconstruct_PT(cfg, imag(PP_SFR{n}), fpsi_bg, maxiter, tol, symmode);
    out.PT_imag{n} = out3.PT;
end
end