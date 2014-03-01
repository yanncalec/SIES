% Lim view identification. Un peu de bidouillage. 

%% Add path
clear all;
close all;
clc;
addpath('~/OOP/');

load ~/Data/smalldico.mat;
% load ~/Data/MSR.mat;

Nv = 128;

%% Data simulation
% load ~/Data/smalldico.mat;

disp('Data simulation...');

% Rigid transform parameters of a standard shape
scl = 1.5; trl = [-0.5,0.5]'; rtn = pi/3;
% scl = 1.75; trl = [-0.5,0.5]'; rtn = pi/3;

mcenter = [0,0]'; 
mradius = 10; 
cfg = acq.Planewave(N0, mcenter, mradius, N0, [1, 2*pi, 2*pi]);

N0 = Nv; % Number of sources (and receivers)

% Scanning frequencies
sclrange = [1, 2];
sfrange = frange./sclrange;

Nsf = floor(Nf/sclrange(2));
sfreq = linspace(sfrange(1), sfrange(2), Nsf);

Data = {};

for n = 1:length(D)
    fprintf('Processing the shape: %s...\n', D{n}.name_str);

    D1 = (D{n}<rtn) * scl + trl;
    
    P = PDE.Helmholtz_R2(D1, cfg, 0, pmtt_bg, pmeb_bg); 
    % figure; plot(P, 'LineWidth', 1); axis image;

    Data{n} = P.data_simulation(sfreq);
    
end

save('~/Data/MSR-fullview-farfield.mat', 'P', 'Data', 'sfreq', 'sfrange', 'Nsf', 'scl', 'trl', 'rtn', '-v7.3');

%% Build dico
Dico_SD = {};
Dico_FFP = {};

% bdwidth = 1; % anglular aperture ~ 3 deg
bdwidth = 3; 

% bdwidth = 5;
% bdwidth = 11;
% bdwidth = 21;

for n=1:length(D)
    [S, G] = dico.Helmholtz.ShapeDescriptorSCT(SCT{n}, Nv, bdwidth);
    Dico_SD{n} = S;
    Dico_FFP{n} = G;
end

%% Plot

% mask = tools.bandiag_mask(Nv, bdwidth);

% idx = 3;
% Sd = squeeze(mean(mean(Dico_SD{idx}, 1),2));

% fig1 = figure; plot(freq, Sd); hold on; grid on;
% xlim([0.4, 4.1]*pi); 

% S = P.MSR2FFP();

% Sm = squeeze(mean(mean(S, 1),2));
% fig1 = figure; plot(sfreq, Sm); hold on; 

% xlim([1.5, 2.1]*pi); grid on;

%% Matching

mask = tools.bandiag_mask(Nv, bdwidth);
% figure; imagesc(mask); axis image;

Nlvl = 5;

Err = zeros(length(D), length(D), Nlvl+1);
Idx = Err;
Scl = Err;

for n=2 %:Nlvl+1
    nlvl = (n-1)/Nlvl
    
    for idx=1:length(D)
        fprintf('Shape %d\n',idx);
        
        Data{idx} = P.add_white_noise(Data{idx}, nlvl);

        A = P.MSR2FFP(Data{idx}.MSR_noisy, sfreq, mask);
        
        S = dico.Helmholtz.ShapeDescriptorFFP(A);

        [t1, t2, t3] = dico.Helmholtz.SCT_matching(S, sfrange, Dico_SD, frange);
        Err(idx,:, n) = t1; 
        Idx(idx,:, n) = t2;
        Scl(idx,:, n) = t3;
    end    
end

%% Plot

N = length(D);
xdata = 1:N;
% xlabels = { 'Ell', 'Flow', 'A', 'Squa',  'E', 'Rect', 'Circ', 'Tria'}; %Dico.name;
xlabels = { 'Ellipse', 'Flower', 'A', 'Square',  'E', 'Rectangle', 'Circle', 'Triangle'}; %Dico.name;

nl=2;
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
set(gca, 'Xtick', xdata, 'XtickLabel', xlabels);

% ylim([0, 140]);

% saveas(fig1, ['../figures/ident_nlvl_0p',num2str(2*(nl-1)),'_limview_wide.eps'],'psc2');
saveas(fig1, ['../figures/ident_nlvl_0p',num2str(2*(nl-1)),'_limview_8deg.eps'],'psc2');

% % nl=3;
Scl1 = squeeze(Scl(:,:,nl));
% fig2 = figure;
% bar(diag(Scl1)-1.5); hold on;
% set(gca, 'Xtick', xdata, 'XtickLabel', xlabels);
% saveas(fig2, ['../figures/sclest_nlvl_0p',num2str(nl-1),'_fullview.eps'],'psc2');
