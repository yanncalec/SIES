%% Add path
clear all;
close all;
clc;
addpath('~/SIES/');

load ~/Data/smalldico.mat;
load ~/Data/MSR.mat;

% ord = (size(Dico.SCT{1},1)-1)/2;

%% Matching 

Nlvl = 5;

Err = zeros(length(D), length(D), Nlvl+1);
Idx = Err;
Scl = Err;

% sfreq0=sfreq;
% sfreq =sfreq0(1:floor(length(sfreq0)/4));
sfreq =sfreq0;

for n=1:Nlvl+1
    nlvl = (n-1)/Nlvl

    for idx=1:length(D)
        S = zeros(Nv, Nv, length(sfreq));
        
        fprintf('Shape %d\n',idx);
        for f=1:length(sfreq)
            fprintf('Freqency %d\n',f);
            Data{idx} = P.add_white_noise(Data{idx}, nlvl);
            MSR = Data{idx}.MSR_noisy{f};
            P.freq = sfreq(f);
            out = P.reconstruct_SCT_analytic(MSR, ord);
            S(:,:,f) = dico.Helmholtz.ShapeDescriptorSCT(out.SCT, Nv);
        end

        [t1, t2, t3] = dico.Helmholtz.SCT_matching(S, [sfreq(1), sfreq(end)], Dico.SD_S, frange);
        Err(idx,:, n) = t1; 
        Idx(idx,:, n) = t2;
        Scl(idx,:, n) = t3;
    end    
end

% fprintf('True shape: %s\n', D{idx}.name_str);
% fprintf('Identified shapes: %s, %s, %s\n', D{idx1(1)}.name_str, D{idx1(2)}.name_str, D{idx1(3)}.name_str);
% fprintf('Errors: %f, %f, %f\n', err1(idx1(1)), err1(idx1(2)), err1(idx1(3)));
% fprintf('Estimations of scaling: %f, %f, %f\n\n', scl1(idx1(1)), scl1(idx1(2)), scl1(idx1(3)));

% Err0 = Err;
% Idx0 = Idx;
% Scl0 = Scl;
%% Plot

N = length(D);
xdata = 1:N;
% xlabels = { 'Ell', 'Flow', 'A', 'Squa',  'E', 'Rect', 'Circ', 'Tria'}; %Dico.name;
xlabels = { 'Ellipse', 'Flower', 'A', 'Square',  'E', 'Rectangle', 'Circle', 'Triangle'}; %Dico.name;

nl=4;
Err1 = squeeze(Err(:,:,nl)); % Each row corresponds to the matching results for one shape
toto = squeeze(Idx(:,:,nl));
Idx1 = zeros(size(Err1));

for n=1:N
    Idx1(n, toto(n,1)) = 1;
end

fig1=figure;
bar(Err1, 'w'); hold on;
bar(eye(size(Err1)).*Err1,'r');

bar(Idx1.*Err1,'g');
ylim([0,750]);
set(gca, 'Xtick', xdata, 'XtickLabel', xlabels);

saveas(fig1, ['../figures/ident_nlvl_0p',num2str(2*(nl-1)),'_fullview.eps'],'psc2');

% nl=3;
Scl1 = squeeze(Scl(:,:,nl));
fig2 = figure;
bar(diag(Scl1)-1.5); hold on;
set(gca, 'Xtick', xdata, 'XtickLabel', xlabels);
saveas(fig2, ['../figures/sclest_nlvl_0p',num2str(2*nl-1)),'_fullview.eps'],'psc2');
