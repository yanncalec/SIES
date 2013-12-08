function Match = identification_xPT_from_multifreq_data(data, freqidx, cfg, ord, maxiter, tol, symmode, Dico, I1, I2, PT, PT_imag)
% Match = identification_xPT_from_multifreq_data(data, freqidx, cfg, ord, maxiter, tol, symmode, Dico, I1, I2, PT, PT_imag)
% Identify a target in a GPT-based dictionary (Shape-Descriptors or PTs).
%
% Inputs:
% data: noisy data, output of the function data_simulation of an ElectricFish object.
% freqidx: the indexes of frequency to be used for matching
% cfg: acquisition configuration
% ord: order of CGPT reconstruction and SD matching
% maxiter: maximum number of iteration for lsqr iteration
% tol: tolerance of lsqr algorithm
% symmode: force the solution to be symmetric
% Dico: a Dictionary object
% I1, I2: SD dictionary obtained from Dico
% PT: PT dictionary obtained from Dico
% PT_imag: imaginary part of PT
%
% Outputs:
% Match: multi-frequency matching

    if (isempty(I1) || isempty(I2)) && isempty(PT) && isempty(PT_imag)
        error('One of the dictionaries (Shape-descriptor, PT or imag(PT)) must be given!');
    end

    % if length(data.freq(freqidx)) == 1
    %     error('Data must be multifrequency!');
    % end
    
    %% Reconstruction of CGPT and PT for all frequencies

    MSR = data.MSR_noisy(freqidx);
    Cur = data.Current_noisy(freqidx);
    PP_SFR = data.PP_SFR_noisy(freqidx);
    
    out = PDE.ElectricFish.reconstruction_xPT_from_multifreq_data(MSR, Cur, PP_SFR, data.fpsi_bg, cfg, ord, maxiter, tol, symmode);    

    %% Multi-frequency matching
    if length(freqidx) > 1
        for n = 1:length(freqidx)
            J1{n} = out.CGPT{n}.I1;
            J2{n} = out.CGPT{n}.I2;
        end

        if ~isempty(I1) && ~isempty(I2)
            Match.SD = Dico.MF_SD_Matching(J1, J2, I1, I2, ord);
        end
        if ~isempty(PT)
            Match.PT = Dico.MF_PT_Matching(out.PT, PT);
        end
        if ~isempty(PT_imag)
            Match.PT_imag = Dico.MF_PT_Matching(out.PT_imag, PT_imag);
        end
    else
        if ~isempty(I1) && ~isempty(I2)
            Match.SD = Dico.SD_Matching(out.CGPT{1}.I1, out.CGPT{1}.I2, I1{1}, I2{1}, ord);
        end
        if ~isempty(PT)
            Match.PT = Dico.PT_Matching(out.PT{1}, PT{1});
        end
        if ~isempty(PT_imag)
            Match.PT_imag = Dico.PT_Matching(out.PT_imag{1}, PT_imag{1});
        end        
    end
end


