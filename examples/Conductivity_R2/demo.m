%% Demo of the Conductivity_R2 class
% This script shows how to use |PDE.Conductivity_R2| class for data simulation
% and reconstruction of GPTs.

%% Add path
clear all;
close all;
clc;
addpath('../../');
%% Definition of small inclusions

%%
% Initialize an object of |C2boundary|

B = shape.Ellipse(1,1/2,2^9);
% B = shape.Flower(1/2, 1/2, 2^10);
% B = shape.Triangle(1/2, pi*0.8, 2^10);
% B = shape.Rectangle(1, 1/2, 2^10);
% B = shape.Banana(5/2, 1, [0,10]', 0, 0, 2^10);
% B = shape.Imgshape('~/Data/images/Letters/R.png', 2^10);

%%
% Make (multiple) inclusion(s)
D{1}=(B<(0.2*pi))*0.5 + 0.25*[1,1]';
% D{2}=B*0.2 + 0.3*[-1,-1]';
cnd = [10, 10]; 
pmtt = [5, 5];

% D{1}=B;
% % D{1}=(B<(0.5*pi))*0.5+[1,1]';
% cnd = [0.5]; 
% pmtt = [0];

%% Set up an environment for experience
% The sources/receptors are distributed on a circle whose center is closed
% to the mass of center of the inclusion, up to a small offset.

%%
% Make the acquisition configuration with the class |acq.Coincided|.
cfg = acq.Coincided([0,0]', 3, 50, [1, 2*pi, 2*pi], 0);
% cfg = acq.Coincided([0,0]', 3, 10, [5, 0.2*pi, 2*pi], 0);

P = PDE.Conductivity_R2(D, cnd, pmtt, cfg); 

fig=figure; plot(P, 'LineWidth', 1); axis image;

%% Simulation of the MSR data
freqlist = linspace(0, 0.1, 3);
tic
data = P.data_simulation(freqlist);
toc

%%
% Calculate and plot potiential fields
sidx = 1;
[F, F_bg, SX, SY, mask] = P.calculate_field(0.01, sidx, [0,0]', 6, 100);
P.plot_field(sidx, F, F_bg, SX, SY, 100);

%% Reconstruction of CGPT
%%
% maximum order of the reconstruction
ord = 4;
symmode = 1;

%%
% Compute first the theoretical value of CGPT
for f=1:length(freqlist)
    lambda = asymp.CGPT.lambda(cnd, pmtt, freqlist(f));
    M{f} = asymp.CGPT.theoretical_CGPT(D, lambda, ord);
end

%%
% add white noise
nlvl = 0.1; 
data = P.add_white_noise(data, nlvl);

%%
% Reconstruct CGPT and show error
fprintf('Norm of the difference between theoretical and reconstructed CGPT matrix at different frequencies:\n');

out = {};
for f=1:length(freqlist)
    MSR = data.MSR{f};
    % out{f} = P.reconstruct_CGPT(MSR, ord, 100000, 1e-10, symmode);
    out{f} = P.reconstruct_CGPT_analytic(MSR, ord);

    % out{f}.res/norm(MSR,'fro')
    norm(M{f} - out{f}.CGPT, 'fro')
end

%%
% The matrix of CGPTs M is symmetric:
fprintf('Difference between theoretical and reconstructed CGPT matrix at the frequency %f:\n', freqlist(2));
M{2}-M{2}.'
