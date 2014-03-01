%% Demo of the Helmholtz class
% This script shows how to use |PDE.Helmhotz_R2| class for data simulation and reconstruction of
% scattering coefficients.

% It's important to know the relations between the frequency, the reconstruction order and the
% truncation error. Higher is the frequency, more is the information in SCT matrix so higher is
% the order. The ellipse is special because among the simple shapes it is the most informative
% one (SCT matrix decays the most slowly). For an ellipse of size 1, with freq=4*pi, it seems
% that the order should be >= 30 to give neglectable (~1e-6) reconstruction error.

%% Add path
clear all;
close all;
clc;
addpath('~/OOP/');

%% Definition of the small inclusion

delta = 1 ; % diameter of the standard shape
%%
% Initialize an object of |C2boundary|

B = shape.Ellipse(delta,delta/2,[0,0]',0,2^10);
% B = shape.Flower(delta/2,delta/2,[0,0]',0,2^10,5,0.4,0.5); 
% B = shape.Triangle(delta/2, pi*0.8, 2^10);
% B = shape.Rectangle(delta,0.5*delta,[0,0]',0,2^10);
% % or load image from file
% B = shape.Imgshape('../images/Letters/R.png', delta, delta, 2^10);
% figure; plot(B); axis image;

%%
% Set the conductivity and the permittivity of the inclusion
B.cnd = 3; B.pmtt = 3; B.pmeb = 3;

%%
% Scaling, translation and rotation parameters:
scl = 1.5; trl = [-0.5,0.5]'; rtn = pi/3;
% scl = 1; trl = [0,0]'; rtn = pi/3;

%%
% The true inclusion _D_ is the shape _B_ after some rigid transform and
% scaling. Apply these transforms:
D = (B < rtn) * scl + trl; % or D = B.t_s_r(trl, scl, rtn);

%% Set up an environment for experience
% The sources/receptors are distributed on a circle whose center is closed
% to the mass of center of the inclusion, up to a small offset.

%%
% offset of the measurement center
% vh = D.diameter/2*(rand(2,1)-0.5)
% vh= [0,0]';
vh = [-0.25,0.25]';

%%
% center of the measurement circle
mcenter = D.center_of_mass + vh; 
%%
% radius of the measurement circle
mradius = 2.5; % D.diameter*2; 

%%
% Make the acquisition configuration. With the class |acq.Planewave|, we make equally distributed
% sources/receivers on a circle of center |mcenter| with radius |mradius|.

N0 = 91; % Number of sources (and receivers)
cfg = acq.Planewave(N0, mcenter, mradius, N0, [1, 2*pi, 2*pi]);

%%
freq = 4*pi; % Working frequency of sources
pmtt_bg = 1; 
pmeb_bg = 1;

P = PDE.Helmholtz_R2(D, cfg, freq, pmtt_bg, pmeb_bg); 

figure; plot(P, 'LineWidth', 1); axis image;

% xs = cfg.src(1); Sx = linspace(-3,3,100); 
% F = PDE.Conductivity_R2.solve_forward(D, 0, xs, Sx, Sx);
% figure; imagesc(F); colorbar()

%% Simulation of the MSR data
% freqlist = linspace(0,1,10);
tic
data = P.data_simulation();
toc

%%
% % Calculate the field and plot it. 
% [F, F_bg, SX, SY, mask] = P.calculate_field(3, cfg.src(1), 2*D.diameter, 64);
% P.plot_field(F, F_bg, SX, SY, 10, 'LineWidth', 2);

%%
% add white noise
nlvl = 0.1; 
data = P.add_white_noise(data, nlvl);
MSR = data.MSR{1}; 
% MSR = data.MSR_noisy{1}; 

%% Compute the SCT
%%
% maximum order of the reconstruction
% ord = floor((N0-1)/2); 
ord = 10;

disp('Computing the SCT matrix of the reference shape...');
% B.cnd = 8; B.pmtt = 8; B.pmeb = 8;
tic
W0 = B.SCT(ord, freq, pmtt_bg, pmeb_bg); 
toc

E = D - cfg.center;
disp('Computing the SCT matrix of the true shape (centered) ...');
tic
W1 = E.SCT(ord, freq, pmtt_bg, pmeb_bg); 
toc

figure; imagesc(abs(W0)); colorbar()
figure; imagesc(abs(W1)); colorbar()

%% Shape descriptors
Nv = 512;
[S0, G0] = dico.Helmholtz.ShapeDescriptorSCT(W0, Nv);
[S1, G1] = dico.Helmholtz.ShapeDescriptorSCT(W1, Nv);

% figure; imagesc(S0); colorbar()
% figure; imagesc(S0-S1); colorbar()
% figure; imagesc(([S1 S1; S1 S1])); colorbar()

disp('Error of the shape descriptor between the true and the reference shape (big if scl<>1):')
norm(S0-S1,'fro')/norm(S0,'fro')

% norm(G0-G1,'fro')/norm(G0,'fro')
% norm(W0-W1,'fro')/norm(W0,'fro')

% figure; imagesc((G0)); colorbar()
% figure; imagesc(G0-flipud(fliplr(G0))); colorbar()

%% Reconstruct SCT and show error
out = P.reconstruct_SCT_analytic(MSR, ord);

% out1 = P.reconstruct_SCT(MSR, ord);
% norm(out0.SCT0-out1.SCT0,'fro')/norm(out1.SCT0)

% figure; imagesc(abs(out.SCT0)); colorbar() % Without post-processing
% figure; imagesc(abs(out.SCT)); colorbar() % With post-processing

disp('Reconstruction error without post-processing:')
norm(out.SCT0-W1,'fro')/norm(W1,'fro')

disp('Reconstruction error with post-processing:')
norm(out.SCT-W1,'fro')/norm(W1,'fro')

disp('Relative truncation error:')
out.res/norm(MSR,'fro')

% W2 = out.SCT;
[S2, G2] = dico.Helmholtz.ShapeDescriptorSCT(out.SCT, Nv);

disp('Error of the shape descriptor between the reconstructed and the true shape:')
norm(S1-S2,'fro')/norm(S1,'fro')

disp('Error of the shape descriptor between the reconstructed and the reference shape (big if scl<>1):')
norm(S0-S2,'fro')/norm(S0,'fro')

