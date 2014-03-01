% Lim view identification. Un peu de bidouillage. 

%% Add path
clear all;
close all;
clc;
addpath('~/OOP/');

load ~/Data/smalldico.mat;
load ~/Data/MSR.mat;

%% Get the far field pattern
Nv = 128;
mask = tools.bandiag_mask(Nv, 20);
DicoF.SD_S = {};

for idx=1:length(D)
    W = Dico.SCT{idx};
    F = zeros(Nv,Nv,Nf);
    
    for f=1:Nf
        F0=dico.Helmholtz.farfieldpattern(W(:,:,f), Nv);
        F(:,:,f) = F0.*mask;
    end

    [toto, ~] = dico.Helmholtz.ShapeDescriptorSCT(F, Nv);
    DicoF.SD_S{idx} = toto;
end

%% Matching
toto1=find(freq>=0.75*pi);
toto2=find(freq<=3*pi);
sfidx = [toto1(1), toto2(end)]

Nlvl = 10;

Err = zeros(length(D), length(D), Nlvl+1);
Idx = Err;
Scl = Err;

for n=1:2
    nlvl = (n-1)/Nlvl

    for idx=1:length(D)
        fprintf('Shape %d\n',idx);
        S = zeros(Nv, Nv, length(sfreq));
        
        for f=sfidx(1):sfidx(2)
            %fprintf('Freqency %d\n',f);
            toto = DicoF.SD_S{idx};
            S(:,:,f-sfidx(1)+1) = tools.add_white_noise(toto(:,:,f), nlvl, 0);
        end

        [t1, t2, t3] = dico.Helmholtz.SCT_matching(S, [freq(sfidx(1)), freq(sfidx(2))], DicoF.SD_S, frange);
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

nl=1;
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

% ylim([0,750]);

saveas(fig1, ['../figures/ident_nlvl_0p',num2str(nl-1),'_limview_srange.eps'],'psc2');

% % nl=3;
Scl1 = squeeze(Scl(:,:,nl));
% fig2 = figure;
% bar(diag(Scl1)-1.5); hold on;
% set(gca, 'Xtick', xdata, 'XtickLabel', xlabels);
% saveas(fig2, ['../figures/sclest_nlvl_0p',num2str(nl-1),'_fullview.eps'],'psc2');
