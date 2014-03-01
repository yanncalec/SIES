%% Plot the acquisition systems and the SVD of the operators

%% Add path
clear all;
% close all;
clc;
addpath('~/OOP/');

%% Definition of the small inclusion

delta = 1 ; % diameter of the standard shape
%%
% Initialize an object of |C2boundary|

% B = shape.Ellipse(delta,delta/2,[0,0]',0,2^10);

B = shape.Flower(delta/2,delta/2,[0,0]',0,2^10,5,0.4,0); 
% B = shape.Triangle(delta/2, pi*0.8, 2^10);
% B = shape.Rectangle(delta,0.5*delta,[0,0]',0,2^10);
% % or load image from file
% B = shape.Imgshape('../images/Letters/R.png', delta, delta, 2^10);
% figure; plot(B); axis image;

%%
% Set the conductivity and the permittivity of the inclusion
B.cnd = 3; B.pmtt = 3; B.pmeb = 3;

%% Set up an environment for experience
% The sources/receptors are distributed on a circle whose center is closed
% to the mass of center of the inclusion, up to a small offset.

%%
% center of the measurement circle
mcenter = [0,0]';
%%
% radius of the measurement circle
mradius = 10;

%%
% Make the acquisition configuration. With the class |acq.Planewave|, we make equally distributed
% sources/receivers on a circle of center |mcenter| with radius |mradius|.
N0 = 512; % Number of sources (and receivers)
cfg = acq.Planewave(N0, mcenter, mradius, N0, [1, 2*pi, 2*pi]);

%%
% Scaling, translation and rotation parameters:
scl = 1; trl = [0,0]'; rtn = 0;
% trl = D.diameter/2*(rand(2,1)-0.5); 

%%
% The true inclusion _D_ is the shape _B_ after some rigid transform and
% scaling. Apply these transforms:
D = (B < rtn) * scl + trl; % or D = B.t_s_r(trl, scl, rtn);

%%
freq = 2*pi; % Working frequency of sources
pmtt_bg = 1; 
pmeb_bg = 1;

P = PDE.Helmholtz_R2(D, cfg, freq, pmtt_bg, pmeb_bg); 
fig = figure; plot(P); axis image; 

tic
data = P.data_simulation();
toc

ord = 30;
W0 = D.SCT(ord, freq, pmtt_bg, pmeb_bg);

%% Far field pattern
Nv = N0;
G = dico.Helmholtz.farfieldpattern(W0, Nv);
G = G / sqrt(8*pi*P.wavenb_bg*mradius);

MSR = data.MSR{1};
V = abs(MSR);

% G1 = fliplr(G0);

% tt=floor(Nv/4)*3;
% G=circshift(G1, [tt, tt]);;

fig1=figure; imagesc(G); axis image; colorbar();

fig2=figure; imagesc(V); axis image; colorbar();

saveas(fig1,'../figures/ffpattern.eps','psc2');
saveas(fig2,'../figures/Vsr.eps','psc2');

fig3=figure; imagesc(abs(G-V)); axis image; colorbar();
norm(G-V,'fro')/norm(V,'fro')