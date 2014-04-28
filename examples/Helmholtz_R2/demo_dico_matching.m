%% Demo of dictionary matching with the Helmholtz class
% This script shows how to use |PDE.Helmhotz_R2| class for dictionary matching

%% Add path

clear all;
close all;
clc;
addpath('../../');

%% Load the dictionary and simulated data. 
% Recall that the inclusion belongs to the dictionary of shapes up to some rigid transform and
% scaling. See the script demo_make_dico.m and demo_data_sim.m

load ~/Data/dico/SCT/smalldico_SCT.mat; % structure Dico
load ~/Data/measurements/SCT/MSR.mat; % Structrure Data

%% Recover some constants 
ord = Dico.ord; % order of the dictionary
sfreq = Data.sfreq; % scanning frequency
Nf = length(Dico.freq); % number of frequency in the dictionary
frange = Dico.frange; % dictionary frequency range
sfrange = Dico.sfrange; % range of the scaning frequency

Nv = Dico.Nv; % resolution of the SCT dictionary
lendico = length(Dico.B); % number of shapes in the dictionary

%%
% Permittivity and permeability
pmtt = Dico.pmtt;
pmeb = Dico.pmeb;
pmtt_bg = Dico.pmtt_bg;
pmeb_bg = Dico.pmeb_bg;

%%
% Initialize an environment for reconstruction of SCT, the inclusion
% Dico.B{1} here can be arbitrary and has no effect.
P = PDE.Helmholtz_R2(Dico.B{1}, pmtt, pmeb, pmtt_bg, pmeb_bg, Data.cfg);

%% Reconstruction of SCT and matching in the dictionary

%%
% level of noise
nlvl = 0.2; 

%%
% For each shape in the dictionary, do the reconstruction and matching. Each row in the matrix Err, Idx, Scl
% contains the result of one shape in the dictionary.

Err = zeros(lendico); % Euclidian distance between one SCT shape descriptor and the dictionary
Idx = zeros(lendico); % Similarity (index of the dictionary elements) sorted in decreasing order
Scl = zeros(lendico); % Estimation of the scaling constant

for n=1:lendico
    S = zeros(Nv, Nv, length(sfreq));
    
    % Iteration on the frequency. The invariant of scaling is possible with
    % multifrequency data.    
    for f=1:length(sfreq)
        % Add noise
        toto = P.add_white_noise(Data.out{n}, nlvl);
        MSR = toto.MSR_noisy{f};

        % Analytical reconstruction
        out = P.reconstruct_SCT_analytic(MSR, sfreq(f), ord);

        % Compute the SCT shape descriptor using reconstructed
        % coefficients. Recall that this descriptor is invariant to
        % translation and rotation only.
        S(:,:,f) = dico.SCT.ShapeDescriptor_SCT(out.SCT, Nv);
    end

    % Dictionary matching
    [t1, t2, t3] = dico.SCT.SCT_matching(S, sfrange, Dico.SD_S, frange);
    
    Err(n,:) = t1;
    Idx(n,:) = t2;
    Scl(n,:) = t3;    
end

%% Interpretation of the result

%% 
% We show in a bar figure the similarity between dictionary shape descriptors and the one reconstructed from data.

fig1= figure; 
bar(Err, 'facecolor', 'none'); hold on; 
set(gca, 'XTickLabel', Dico.names, 'XTick',1:lendico);

%%
% Identification result is plotted in green, if the result is wrong, the
% true shape is marked in red.
bar(eye(lendico).*Err, 'r'); 

toto=eye(lendico); 
idx = Idx(:,1); 
bar(toto(idx, :).*Err, 'g'); 

%%
% Estimation of scaling for each shape in the dictionary. The bar is the error in the estimation.

disp('True scaling constant:');
Data.scl

fig2 = figure;
bar(diag(Scl)-Data.scl); hold on;
set(gca, 'XTickLabel', Dico.names, 'XTick',1:lendico);


