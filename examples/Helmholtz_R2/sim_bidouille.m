%% Bidouille, don't run

load /media/620b85d2-25b8-4ab1-b5c2-693a81acafae/Data/dico/Helmholtz/smalldico.mat;

% Range of scaling
sclrange = [1/2, 2];

% Range of scanning frequencies
sfrange = [2, 10]*pi;

frange = sclrange .* sfrange;

% Frequency range of the dictionary
Nf = 256; 
freq = linspace(frange(1), frange(2), Nf);

ord = 50;
Dico1 = {};

Nv = 256; % dimension of the shape descriptor

for m=1:length(D)
    fprintf('Processing the shape %d...\n', m);

    W = zeros(2*ord+1, 2*ord+1, Nf);
    S = zeros(Nv, Nv, Nf);
    G = zeros(Nv, Nv, Nf);
    
    for n=1:Nf
        W(:,:,n) = Dico{m}.SCT{n};
        S(:,:,n) = Dico{m}.SD_S{n};
        G(:,:,n) = Dico{m}.SD_G{n};
    end

    Dico1{m}.SCT = W; % SCT
    Dico1{m}.SD_S = S; % Shape descriptor S
    Dico1{m}.SD_G = G; % Shape descriptor G    
end

Dico = Dico1;
clear ans fmax fmin fname scanfreq_range scl_range W S G m n Dico1 freqlist;

save('~/Data/smalldico.mat', '-v7.3');

%% Bidouille 2
% load /media/620b85d2-25b8-4ab1-b5c2-693a81acafae/Data/dico/Helmholtz/MSR.mat;

disp('Data simulation...');

N0 = 181; % Number of sources (and receivers)

% Rigided transform parameters of a standard shape
scl = 1.5; trl = [-0.25,0.25]'; rtn = pi/3;

% Scanning frequencies
sfreq = linspace(sfrange(1), sfrange(2), 64);

%%
% center of the measurement circle
mcenter = [0,0]'; 

%%
% radius of the measurement circle
mradius = 2.5; 

%%
% Make the acquisition configuration. With the class |acq.Planewave|, we make equally distributed
% sources/receivers on a circle of center |mcenter| with radius |mradius|.
cfg = acq.Planewave(N0, mcenter, mradius, N0, [1, 2*pi, 2*pi]);

for n = 3
    fprintf('Processing the shape: %s...\n', D{n}.name_str);

    D1 = (D{n}<rtn) * scl + trl;
    
    P = PDE.Helmholtz_R2(D1, cfg, 0, pmtt_bg, pmeb_bg); 
    figure; plot(P, 'LineWidth', 1); axis image;

    Data = P.data_simulation(sfreq);
end

save('~/Data/MSR_D3.mat', 'P', 'Data', 'sfreq', 'scl', 'trl', 'rtn', '-v7.3');


%% Compute the SD of a transformed shape
clc;
clear all;
close all;
addpath('~/OOP/');
imagepath = '~/Data/images/Letters/';

load '~/Data/smalldico.mat';
load '~/Data/MSR_D3.mat';

WD = P.D{1}.SCT(ord, sfreq, pmtt_bg, pmeb_bg);
[SD, GD] = dico.ShapeDescriptorSCT(WD, Nf);

save('~/Data/WD.mat', 'WD','SD','GD','-v7.3');

