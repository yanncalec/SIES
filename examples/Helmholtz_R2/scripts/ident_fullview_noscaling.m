%% Identification without scaling by comparing directly the shape descriptor
% I tried to find the resemblance between shapes but visiblely 
% Ellipse ~ circle
% Flower ~ circle
% A ~ E
% Square ~ circle
% E ~ A
% Rectangle ~ square
% Circle ~ flower
% Triangle ~ E

% Good news is the identification is very robust to noise.

%% Add path
clear all;
close all;
clc;
addpath('~/OOP/');

load ~/Data/smalldico.mat;

fidx = 1;
Nv = 512;

for n=1:length(D)
    toto = SCT{n}; 
    [S, G] = dico.Helmholtz.ShapeDescriptorSCT(toto(:,:,fidx), Nv);
    Dico_SD{n} = S;
    Dico_FFP{n} = G;
end

%% Data simulation
% load ~/Data/smalldico.mat;

disp('Data simulation...');

sfreq = freq(fidx);

N0 = 91; % Number of sources (and receivers)

Data = {};

scl = 1; trl = [-0.5,0.5]'; rtn = pi/3;
mcenter = [0,0]';  mradius = 3; 

cfg = acq.Planewave(N0, mcenter, mradius, N0, [1, 2*pi, 2*pi]);

for n = 1:length(D)
    fprintf('Processing the shape: %s...\n', D{n}.name_str);

    D1 = (D{n}<rtn) * scl + trl;
    
    P = PDE.Helmholtz_R2(D1, cfg, sfreq, pmtt_bg, pmeb_bg); 
    figure; plot(P, 'LineWidth', 1); axis image;

    Data{n} = P.data_simulation(sfreq);
end

%% Matching 

ord = 30;
Nlvl = 5;

Err = zeros(length(D), length(D), Nlvl+1);
Idx = Err;

for n=1:Nlvl+1
    nlvl = (n-1)/Nlvl*2;

    for idx=1:length(D)
        fprintf('Shape %d\n',idx);

        Data{idx} = P.add_white_noise(Data{idx}, nlvl);
        MSR = Data{idx}.MSR_noisy{1};
        out = P.reconstruct_SCT_analytic(MSR, ord);
        S = dico.Helmholtz.ShapeDescriptorSCT(out.SCT, Nv);

        [t1, t2] = dico.Helmholtz.SCT_matching_noscaling(S, Dico_SD);
        Err(idx,:, n) = t1; 
        Idx(idx,:, n) = t2;
    end    
end

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

set(gca, 'Xtick', xdata, 'XtickLabel', xlabels);

%saveas(fig1, ['../figures/ident_nlvl_0p',num2str((nl-1)),'_fullview_noscl.eps'],'psc2');
