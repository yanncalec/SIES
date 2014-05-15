%% Demo of the Helmholtz_R2 class
% This script shows how to use |PDE.Helmhotz_R2| class for data simulation and reconstruction of
% scattering coefficients. We calculate also the shape descriptor which is
% invariant to the translation and rotation (but not to scaling).

%% Add path
clear all;
close all;
clc;
addpath('../../');

%% Definition of small inclusions

%%
% Initialize an object of |C2boundary|

% B = shape.Ellipse(1,1/2,2^9);
B = shape.Flower(1/2, 1/2, 2^10);
% B = shape.Triangle(1/2, pi*0.8, 2^10);
% B = shape.Rectangle(1, 1/2, 2^10);
% B = shape.Banana(5/2, 1, [0,10]', 0, 0, 2^10);
% B = shape.Imgshape('../images/Letters/R.png', 2^10);

%%
% Make (multiple) inclusion(s). The inclusion D is obtained by rotating and
% translating the reference shape B.
D{1} = (B<(0.2*pi)) + 0.25*[1,1]';

%%
% Set the conductivity and the permittivity of the inclusion and the
% background
pmtt_bg = 1; pmeb_bg = 1;
pmeb = 3; pmtt = 3;

%% Set up an environment for experience
% The sources/receptors are distributed on a circle whose center is closed
% to the mass of center of the inclusion, up to a small offset.

%%
% Make the acquisition configuration. With the class |acq.Planewave|, we make equally distributed
% sources/receivers on a circle.
N0 = 51; % Number of sources (and receivers)
cfg = acq.Planewave([0,0]', 2.5, N0, N0, [1, 2*pi, 2*pi]);

%%
% Initialize an |Helmholtz_R2| object 
P = PDE.Helmholtz_R2(D, pmtt, pmeb, pmtt_bg, pmeb_bg, cfg); 
figure; plot(P, 'LineWidth', 1); axis image;

%% Simulation of the MSR data
freq = 2*pi; % Working frequency of sources

tic
data = P.data_simulation(freq);
toc

%%
% Calculate and plot potiential fields
sidx = 1;
[F, F_bg, SX, SY, mask] = P.calculate_field(0.01, sidx, [0,0]', 6, 100);
P.plot_field(sidx, F, F_bg, SX, SY, 50);

%% Compute the theoretical value of CGPT

%%
% maximum order of the reconstruction
ord = floor((N0-1)/2); 
% ord = 4;

%%
% SCT of the shape |D| and the reference shape |B|
WD = asymp.SCT.theoretical_SCT(D, pmeb, pmtt, pmeb_bg, pmtt_bg, ord, freq);
WB = asymp.SCT.theoretical_SCT(B, pmeb, pmtt, pmeb_bg, pmtt_bg, ord, freq);

%%
% Far field pattern
FP = dico.SCT.farfieldpattern(WB, 512);
fig=figure;
imagesc(abs(FP)); axis image; title('Far filed pattern');
% saveas(fig, '~/farfieldpattern.eps','psc2');

mask = tools.bandiag_mask(512,100);
fig=figure;
imagesc(abs(FP).*mask); axis image; title('Far filed pattern with limited angle of view'); 
% saveas(fig, '~/farfieldpattern_limview.eps','psc2');

%% Reconstruct SCT and show error

%%
% add white noise
nlvl = 0.1; 
data = P.add_white_noise(data, nlvl);

% MSR = data.MSR{1};
MSR = data.MSR_noisy{1};

%%
% Reconstruction by analytical inversion. This works only for equally
% distributed sources and receivers.
out = P.reconstruct_SCT_analytic(MSR, freq, ord);

%%
% Show error
disp('Relative error of reconstructed scattering coefficients:')
norm(out.SCT-WD,'fro')/norm(WD,'fro')

disp('Relative error of data-fitting (least square error):')
out.res/norm(MSR,'fro')

%% Compute the shape descriptor from SCT

%%
% Resolution in frequency domain of the shape descriptor S
Nv = 256; 
%%
% Shape descriptor of the reference shape B
[SB, GB] = dico.SCT.ShapeDescriptor_SCT(WB, Nv);

%%
% Shape descriptor of the reconstruction
[SR, GR] = dico.SCT.ShapeDescriptor_SCT(out.SCT, Nv);

%%
% Shape descriptor of the true shape D
[SD, GD] = dico.SCT.ShapeDescriptor_SCT(WD, Nv);

%%
% Error of SD between the true value and the reconstruction should be small
disp('Error of the shape descriptor between the reconstructed and the true one:')
norm(SR-SD,'fro')/norm(SD,'fro')

%%
% D and B have the same shape descriptor
disp('Error of the shape descriptor between the true and the reference one (big if scl<>1):')
norm(SD-SB,'fro')/norm(SB,'fro')

