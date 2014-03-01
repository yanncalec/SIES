%% Demo of the Helmholtz class
% This script shows how to use |PDE.Helmhotz_R2| class for data simulation and reconstruction of
% scattering coefficients. We calculate also the shape decsriptor which is
% invariant to the translation ans rotation (but not to scaling).

%% Add path
clear all;
close all;
clc;
addpath('~/SIES/');

%% Definition of small inclusions

%%
% Initialize an object of |C2boundary|

B = shape.Ellipse(1,1/2,2^9);
% B = shape.Flower(1/2, 1/2, 2^10);
% B = shape.Triangle(1/2, pi*0.8, 2^10);
% B = shape.Rectangle(1, 1/2, 2^10);
% B = shape.Banana(5/2, 1, [0,10]', 0, 0, 2^10);
% B = shape.Imgshape('../images/Letters/R.png', 2^10);

%%
% Make (multiple) inclusion(s)
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
% sources/receivers on a circle of center |mcenter| with radius |mradius|.
N0 = 51; % Number of sources (and receivers)
cfg = acq.Planewave([0,0]', 2.5, N0, N0, [1, 2*pi, 2*pi]);

%%
freq = pi; % Working frequency of sources

P = PDE.Helmholtz_R2(D, pmtt, pmeb, pmtt_bg, pmeb_bg, cfg); 

figure; plot(P, 'LineWidth', 1); axis image;

% xs = cfg.src(1); Sx = linspace(-3,3,100); 
% F = PDE.Conductivity_R2.solve_forward(D, 0, xs, Sx, Sx);
% figure; imagesc(F); colorbar()

%% Simulation of the MSR data
freq = pi;
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

WD = asymp.SCT.theoretical_SCT(D, pmeb, pmtt, pmeb_bg, pmtt_bg, ord, freq);
WB = asymp.SCT.theoretical_SCT(B, pmeb, pmtt, pmeb_bg, pmtt_bg, ord, freq);

%% Reconstruct SCT and show error

%%
% add white noise
nlvl = 0.1; 
data = P.add_white_noise(data, nlvl);

MSR = data.MSR{1};
% MSR = data.MSR_noisy{1};

out = P.reconstruct_SCT_analytic(MSR, freq, ord);

norm(out.SCT-WD,'fro')/norm(WD,'fro')

disp('Relative truncation error:')
out.res/norm(MSR,'fro')

%% Compute the shape descriptor from SCT

Nv = 256; % Resolution in frequency domain of the shape descriptor S
%%
% Shape descriptor of the reference shape B
[S0, G0] = dico.SCT.ShapeDescriptor_SCT(WB, Nv);

%%
% Shape descriptor of the reconstruction
[S1, G1] = dico.SCT.ShapeDescriptor_SCT(out.SCT, Nv);

%%
% Shape descriptor of the true shape D
[S2, G2] = dico.SCT.ShapeDescriptor_SCT(WD, Nv);

%%
% Error of SD between the true value and the reconstruction should be small
disp('Error of the shape descriptor between the reconstructed and the true shape:')
norm(S1-S2,'fro')/norm(S1,'fro')

%%
% D and B have the same shape descriptor
disp('Error of the shape descriptor between the reconstructed and the reference shape (big if scl<>1):')
norm(S0-S2,'fro')/norm(S0,'fro')

