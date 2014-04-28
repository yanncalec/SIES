%% Demo of the PulseImaging_R2 class
% This script shows how to use |PDE.PulseImaging| class for data simulation
% and reconstruction of time dependent CGPTs.

%% Add path
clear all;
close all;
clc;
addpath('../../');

%% Definition of the small inclusion

delta = 1 ; % diameter of the standard shape
%%
% Make (multiple) inclusions

nbPoints = 2^9; % number of discritization points on the boundary

% Initialize an object of |shape.C2boundary|
B = shape.Flower(delta/2, delta/2, nbPoints, 5, 0.4, 0);
D{1}=(B<(0.3*pi))*0.2+0.2*[1,1]'; % rotation, scaling and translation to get the first inclusion

cnd = [10]; % conductivity values of two inclusions
pmtt = [5]; % permittivity values

% One can also add multiple inclusions
% B = shape.Ellipse(delta,delta/2,nbPoints);
% D{2}=(B<(0.3*pi))*0.2+0.2*[-1,-1]'; % second inclusion

% cnd = [10, 5]; % conductivity values of two inclusions
% pmtt = [0.1, 0.2]; % permittivity values

%% Set up an environment for experience
% The sources/receptors are distributed on a circle whose center is closed
% to the mass of center of the inclusion, up to a small offset.

cfg = acq.Coincided([0,0]', 1, 10, [1, 2*pi, 2*pi], 0, [1,-1]);

[waveform,dt] = PDE.PulseImaging_R2.make_pulse(100);
% figure; plot(waveform);

%%
% Initialize an environment by passing the fish, the inclusions, the
% configuration, and the physical constants etc.

P = PDE.PulseImaging_R2(D, cnd, pmtt, waveform, dt, cfg);
% figure; plot(P, 'LineWidth', 1); axis image;

%% Simulation of data
% The output is a structure containing many fields, see the function
% data_simulation().

Ntime = 50; % time interval length

tic
data = P.data_simulation(Ntime);
toc

%% Plot the potential fields
% *Bugs*: some grid points are on the fish's skin which may cause numerical
% instability. 
%%
% Calculate the field and plot it. 
sidx = 1; % source index to be calculated
[F, F_bg, SX, SY] = P.calculate_field(100, sidx, [0,0]', 3, 100);

% P.plot_field(sidx, F, F_bg, SX, SY, 100, 0, '-g','LineWidth', 1.1);
close all;
figure;
for t=1:100
    imagesc(F{t}-F_bg{t}); axis image; colorbar(); 
end
