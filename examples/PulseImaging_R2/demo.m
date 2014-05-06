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

nbPoints = 2^10; % number of discritization points on the boundary

% Initialize an object of |shape.C2boundary|
B = shape.Flower(delta/2, delta/2, nbPoints, 5, 0.4, 0);
% D{1}=(B<(0.3*pi))*0.2+0.2*[1,1]'; % rotation, scaling and translation to get the first inclusion
D{1}=(B<(0.3*pi))*0.4; % rotation, scaling and translation to get the first inclusion

cnd = [3]; % conductivity values of two inclusions
pmtt = [0.5]; % permittivity values

% One can also add multiple inclusions
% B = shape.Ellipse(delta,delta/2,nbPoints);
% D{2}=(B<(0.3*pi))*0.2+0.2*[-1,-1]'; % second inclusion

% cnd = [10, 5]; % conductivity values of two inclusions
% pmtt = [1, 2]; % permittivity values

%% Set up an environment for experience
% The sources/receptors are distributed on a circle whose center is closed
% to the mass of center of the inclusion, up to a small offset.

cfg = acq.Coincided([0,0]', 1, 10, [1, 2*pi, 2*pi], 0, [1,-1]);

Ntime = 1000; % time interval length
[waveform,dt] = PDE.PulseImaging_R2.make_pulse(3,Ntime);
figure; plot(waveform);

%%
% Initialize an environment by passing the fish, the inclusions, the
% configuration, and the physical constants etc.

P = PDE.PulseImaging_R2(D, cnd, pmtt, waveform, dt, cfg);
figure; plot(P, 'LineWidth', 1); axis image;

%% Simulation of data
% The output is a structure containing many fields, see the function
% data_simulation().

tic
data = P.data_simulation(Ntime);
toc

%% Plot the potential fields and make a movie
%%
% Calculate the field and plot it. 
sidx = 1; % source index to be calculated
[F, F_bg, SX, SY] = P.calculate_field(Ntime, sidx, [0,0]', 3, 200);

% P.plot_field(sidx, F{20}, F_bg{20}, SX, SY, 20); % , '-g','LineWidth', 1.5);

%%
% Make a movie

vmin=0; vmax=0;
for t=1:Ntime
    Img{t} = F{t}-F_bg{t}; %perturbation
    vmin = min(vmin, min(min(Img{t})));
    vmax = max(vmax, max(max(Img{t})));
end

nbLine = 30;
colormap(jet);
fig=figure;

for t=1:Ntime
    % fig = figure; 
    % contourf(F{t}-F_bg{t}, nbLine); axis image; colorbar(); 
    imagesc(F{t}-F_bg{t}); axis image; colorbar(); 
    caxis manual;    
    caxis([vmin, vmax]);
    saveas(fig,['img/',num2str(t),'.png']);
end

vidObj = VideoWriter('img/movie.avi');
vidObj.FrameRate=10;
open(vidObj);
for t=1:Ntime
    img = imread(['img/',num2str(t),'.png']);
    writeVideo(vidObj, img);
end
close(vidObj);
