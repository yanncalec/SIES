function out = identification_xPT_one_target_many_exps(data0, freqidx, cfg, nlvl, Nexp, Dico, I1, I2, PT, PT_imag, ord)
% out = identification_xPT_one_target_many_exps(data0, freqidx, cfg, nlvl, Nexp, Dico, I1, I2, PT, PT_imag, ord) 
%
% This function realizes Nexp different noise samples of a same level, and
% for each realization it calls identification_xPT_from_multifreq_data. It
% should be used at the top-level of a numerical experiment.
%
% Inputs:
% data0: noiseless (mutli) frequency data corresponding to one target
% cfg:
% nlvl: noise level
% Nexp: number of noise realization
% Dico: object of dico.Dictionary
% I1, I2: dictionary of Shape Descriptors
% PT, PT_imag: dictionary of PT and imag(PT)
% ord: order for SD matching

%% Reconstruction and identification with multiple realizations of noise

maxiter = 5000;
tol = 1e-5;
symmode = 1; % force the solution to be symmetric 

% Multi-frequency matching error
out.SD.err = zeros(Nexp, Dico.size);
out.PT.err = zeros(Nexp, Dico.size);
out.PT_imag.err = zeros(Nexp, Dico.size);

for n = 1:Nexp
    % add white noise to data, one realization only
    data = PDE.ElectricFish.add_white_noise(data0, nlvl);

    % reconstruction and identification
    Match{n} = PDE.ElectricFish.identification_xPT_from_multifreq_data(data, freqidx, cfg, ord, maxiter, tol, symmode, Dico, I1, I2, PT, PT_imag);
    % [SMatch{n}, MMatch{n}] = PDE.ElectricFish.identification_xPT_from_multifreq_data(data, cfg, ord, maxiter, tol, symmode, ...
    %                                                   Dico, I1, I2, PT, PT_imag);

    out.SD.err(n, :) = Match{n}.SD.err;
    out.PT.err(n, :) = Match{n}.PT.err;
    out.PT_imag.err(n, :) = Match{n}.PT_imag.err;
end

% Output

if ~isempty(I1) && ~isempty(I2)
    [~, out.SD.idx] = sort(mean(out.SD.err, 1));
    out.SD.name = Dico.name(out.SD.idx);
end

if ~isempty(PT)
    [~, out.PT.idx] = sort(mean(out.PT.err, 1));
    out.PT.name = Dico.name(out.PT.idx);
end

if ~isempty(PT_imag)
    [~, out.PT_imag.idx] = sort(mean(out.PT_imag.err, 1));
    out.PT_imag.name = Dico.name(out.PT_imag.idx);
end

