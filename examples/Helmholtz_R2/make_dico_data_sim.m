%% Make a small dictionary and simulate data for all shapes in the dictionary and for multiple frequencies.

%%
clc;
clear all;
close all;
addpath('~/OOP/');
imagepath = '~/Data/images/Letters/';

%% Construction of the dictionary
disp('Construction of the dictionary...');

nbPoints = 2^10;
delta = 1;

D={};

% D{1} = shape.Ellipse(delta/2,delta/2,[0,0]',0,nbPoints); % disk
% D{2} = shape.Triangle(delta, pi/3, nbPoints);
% D{3} = shape.Rectangle(delta, delta, nbPoints); % square
% D{4} = shape.Ellipse(delta*1,delta/2,[0,0]',0,nbPoints); % ellipse
% D{5} = shape.Flower(delta/2,delta/2,[0,0]',0,nbPoints,5,0.4,0);
% D{6} = shape.Imgshape([imagepath,'A.png'], delta, delta, nbPoints);
% D{7} = shape.Imgshape([imagepath,'E.png'], delta, delta, nbPoints);
% D{8} = shape.Rectangle(delta, delta/2, nbPoints);

D{1} = shape.Ellipse(delta*1,delta/2,[0,0]',0,nbPoints); % ellipse
D{2} = shape.Flower(delta/2,delta/2,[0,0]',0,nbPoints,5,0.4,0);
D{3} = shape.Imgshape([imagepath,'A.png'], delta, delta, nbPoints);
D{4} = shape.Rectangle(delta, delta, nbPoints); % square
D{5} = shape.Imgshape([imagepath,'E.png'], delta, delta, nbPoints);
D{6} = shape.Rectangle(delta, delta/2, nbPoints);
D{7} = shape.Ellipse(delta/2,delta/2,[0,0]',0,nbPoints); % disk
D{8} = shape.Triangle(delta, pi/3, nbPoints);

for n=1:length(D)
    D{n}.cnd = 3;    
    D{n}.pmtt = 3;
    D{n}.pmeb = 3;
end

% Background
pmtt_bg = 1; 
pmeb_bg = 1;

%% Compute the SCT of a reference shape

% Range of scaling
sclrange = [1/2, 2];

% Range of scanning frequencies
sfrange = [1, 2]*pi;

frange = sclrange .* sfrange;

% Frequency range of the dictionary
Nf = floor((frange(2)-frange(1))/0.05); 
freq = linspace(frange(1), frange(2), Nf);

ord = 30;

Dico.SCT = {};
% Dico.D = D;
% Dico.SD_S = {};
% Dico.SD_G = {};

Nv = 512; % dimension of the shape descriptor

for n=1:length(D)
    fprintf('Processing the shape %d...\n', n);

    W = D{n}.SCT(ord, freq, pmtt_bg, pmeb_bg);
    [S, G] = dico.Helmholtz.ShapeDescriptorSCT(W, Nv);

    Dico.SCT{n} = W; % SCT
    % Dico.SD_S{n} = S; % Shape descriptor S
    % Dico.SD_G{n} = G; % Shape descriptor G    
end

clear W S G n;
save('~/Data/smalldico.mat', '-v7.3');

%% Data simulation
% load ~/Data/smalldico.mat;

disp('Data simulation...');

N0 = 91; % Number of sources (and receivers)

% Rigid transform parameters of a standard shape
scl = 1.5; trl = [-0.5,0.5]'; rtn = pi/3;
% scl = 1.75; trl = [-0.5,0.5]'; rtn = pi/3;

% Scanning frequencies
Nsf = floor(Nf/sclrange(2));
sfreq = linspace(sfrange(1), sfrange(2), Nsf);

Data = {};

%%
% center of the measurement circle
mcenter = [0,0]'; 

%%
% radius of the measurement circle
% mradius = 2.5; 
mradius = 3.5; 

%%
% Make the acquisition configuration. With the class |acq.Planewave|, we make equally distributed
% sources/receivers on a circle of center |mcenter| with radius |mradius|.
cfg = acq.Planewave(N0, mcenter, mradius, N0, [1, 2*pi, 2*pi]);

% Dico1.D = {};
% Dico1.SCT = {};
% Dico1.SD_S = {};
% Dico1.SD_G = {};

for n = 1:length(D)
    fprintf('Processing the shape: %s...\n', D{n}.name_str);

    D1 = (D{n}<rtn) * scl + trl;
    
    P = PDE.Helmholtz_R2(D1, cfg, 0, pmtt_bg, pmeb_bg); 
    % figure; plot(P, 'LineWidth', 1); axis image;

    Data{n} = P.data_simulation(sfreq);
    
    % W = D1.SCT(ord, sfreq, pmtt_bg, pmeb_bg);
    % [S, G] = dico.Helmholtz.ShapeDescriptorSCT(W, Nv);

    % Dico1.D{n} = D1;
    % Dico1.SCT{n} = W; % SCT
    % Dico1.SD_S{n} = S; % Shape descriptor S
    % Dico1.SD_G{n} = G; % Shape descriptor G    
end

save('~/Data/MSR.mat', 'P', 'Data', 'sfreq', 'scl', 'trl', 'rtn', '-v7.3');
