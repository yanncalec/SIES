%% Demo of the Conductivity_R2 class
% This script shows the procedure of dictionary matching using (CGPT) shape descriptors with frequency=0. 
% Remark that if the size of the inclusion is too large then the reconstruction of high order CGPT might be inaccurate due to the large truncation error. On the contrary, if the size is too small then different objects might be ressemble.
% The efficiency of the identification depends on the size of the target
% (or the distance to the transmitters). In a noisy environment (10% of
% noise), one should reduce the field distance and use low order (2 or 3)
% shape descriptors for identification.

%% Add path
clear all;
close all;
clc;
addpath('../../');

%% Load the dictionary
load('~/Data/dico/CGPT/smalldico_CGPT.mat');
mydico = Dico;
lendico = length(mydico.B);

%%
% Names of dictionary elements
names = {};
for n=1:lendico
    names{n} = mydico.B{n}.name_str;
end

%% Simulation and identification for a fixed dictionary element

%%
% Choose one element
B = mydico.B{4};
% figure; plot(B, 'LineWidth', 2); axis image;

%%
% Make inclusion
D{1}=(B<(0.2*pi))*0.75 + 0.25*[1,1]';
cnd = mydico.cnd;
pmtt = mydico.pmtt;

%% Set up an environment for experience
% The sources/receptors are distributed on a circle whose center is closed
% to the mass of center of the inclusion, up to a small offset.

%%
% Make the acquisition configuration with the class |acq.Coincided|.
cfg = acq.Coincided([0,0]', 1.5, 100, [1, 2*pi, 2*pi], 0);
% cfg = acq.Coincided([0,0]', 1.5, 100, [5, 0.2*pi, 2*pi], 0);

P = PDE.Conductivity_R2(D, cnd, pmtt, cfg); 

fig=figure; plot(P, 'LineWidth', 1); axis image;

%% Simulation of the MSR data
freqlist = 0;
tic
data = P.data_simulation(freqlist);
toc

%% Reconstruction of CGPT from the noisy MSR data

ord = 5; % maximum order of the reconstruction
cord = 5; % order of comparison

%%
% add white noise
nlvl = 0.05; 
data = P.add_white_noise(data, nlvl);

%%
% Reconstruct CGPT and show error
out = {};
MSR = data.MSR_noisy{1};
% out{f} = P.reconstruct_CGPT(MSR, ord, 100000, 1e-10, symmode);
out{1} = P.reconstruct_CGPT_analytic(MSR, ord);

%% Dictionary matching
[I1, I2, ~] = dico.CGPT.ShapeDescriptor_CGPT(out{1}.CGPT);
[err, idx] = dico.CGPT.SD_Matching(I1, I2, mydico.I1, mydico.I2, cord);

%%
% Interpretation of the result. We show in a bar figure the similarity
% between dictionary shape descriptors and the one reconstructed from data.

figure;
bar(err{3}, 'b');
set(gca, 'XTickLabel',names, 'XTick',1:lendico)

% for n=1:cord    
%     fprintf('Identification with first %d orders CGPT shape descriptors:\n', n);
%     disp(names(idx{n}));
%     disp('Errors in increasing order:');
%     disp(err{n}');
%     
% end

%% Simulation and identification for the whole dictionary
Err={}; Idx={};

for n=1:lendico
    B = mydico.B{n};
    
    D{1}=(B<(0.2*pi))*0.75 + 0.25*[1,1]';
    
    P = PDE.Conductivity_R2(D, cnd, pmtt, cfg);
    
    data = P.data_simulation(freqlist);
    
    data = P.add_white_noise(data, nlvl);
    
    out = {};
    MSR = data.MSR_noisy{1};
    out{1} = P.reconstruct_CGPT_analytic(MSR, ord);
    
    [I1, I2, ~] = dico.CGPT.ShapeDescriptor_CGPT(out{1}.CGPT);
    [toto1, toto2] = dico.CGPT.SD_Matching(I1, I2, mydico.I1, mydico.I2, cord);
    
    Err{n} = cell2mat(toto1);
    Idx{n} = cell2mat(toto2);
end

%%
% Reshape into 3D array. first dim: index of dictionary element, second dim: order of
% comparison, third dim: group of experiment.

Err3D = reshape(cell2mat(Err), lendico, cord, lendico);
Idx3D = reshape(cell2mat(Idx), lendico, cord, lendico);

%% Plot
err = squeeze(Err3D(:, 2, :)); err = err'; % Take the transpose
idx = squeeze(Idx3D(1, 2, :)); idx = idx';

figure; 
bar(err, 'facecolor', 'none'); hold on; 
set(gca, 'XTickLabel',names, 'XTick',1:lendico);

% Identification result is plotted in green, if the result is wrong, the
% true shape is marked in red.
bar(eye(lendico).*err, 'r'); 
toto=eye(lendico); 
bar(toto(idx,:).*err, 'g'); 



