%% Plot the curve S_D and S_B in function of frequency

%% Add path
clear all;
close all;
clc;
addpath('~/OOP/');

load ~/Data/smalldico.mat;
load ~/Data/MSR.mat;

idx = 2;
S = zeros(Nv, Nv, length(sfreq));

for n=1:length(sfreq)
    MSR = Data{idx}.MSR{n};
    P.freq = sfreq(n);
    out = P.reconstruct_SCT_analytic(MSR, ord);
    S(:,:,n) = dico.Helmholtz.ShapeDescriptorSCT(out.SCT, Nv);
end

%% plot
Sd = squeeze(mean(mean(Dico.SD_S{2}, 1),2));
Sm = squeeze(mean(mean(S, 1),2));

fig1 = figure; plot(freq, Sd); hold on; 
xlim([0.4, 4.1]*pi); grid on;

fidx = find((freq>=1.5*pi).*(freq<=3*pi));
plot(freq(fidx), Sd(fidx), 'r','LineWidth',2);

saveas(fig1,'../figures/SB.eps', 'psc2');

fig2 = figure; plot(sfreq, Sm); grid on;
saveas(fig2,'../figures/SD.eps', 'psc2');
