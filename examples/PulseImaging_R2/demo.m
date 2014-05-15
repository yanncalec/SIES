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
% D{2}=(B<(0.3*pi))*0.2+2*[-1,-1]'; % second inclusion

% cnd = [10, 5]; % conductivity values of two inclusions
% pmtt = [1, 2]; % permittivity values

%% Set up an environment for experience
% The sources/receptors are distributed on a circle whose center is closed
% to the mass of center of the inclusion, up to a small offset.

cfg = acq.Coincided([0,0]', 1, 20, [1, 2*pi, 2*pi], 0, [1,-1]);

%% Show the pulse waveform h and its Fourier transfom H
Tmax = 5; Ntime = 2^10; % time interval length
% [waveform, dt, hfunc, freqform, df, Hfunc] = PDE.PulseImaging_R2.make_pulse(2,Ntime);
[waveform, dt, freqform, df] = PDE.PulseImaging_R2.make_pulse(Tmax, Ntime);
figure; plot(linspace(0, Tmax, Ntime), waveform); title('h'); xlabel('time');

Nfreq = length(freqform);
Fmax = (Nfreq-1) * df;
figure; plot(linspace(0, Fmax, Nfreq), abs(freqform)); title('H'); xlabel('frequency');

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

%% Reconstruction of time-dependent CGPT
%%
% add white noise
nlvl = 0.1; 
data = P.add_white_noise(data, nlvl);

%%
% maximum order of the reconstruction
ord = 1;
symmode = 1;

%%
% Compute theoretical values of time-dependent CGPTs
[CGPTt, dt0] = asymp.CGPT.theoretical_CGPT_time(D, cnd, pmtt, ord, freqform, df);

CGPTt0 = asymp.CGPT.CGPT_time_truncation(CGPTt, dt0, Tmax);

Ntime0 = size(CGPTt0,3);
nrm0 = zeros(1,Ntime0);

for t=1:Ntime0
    %    nrm0(t) = norm(CGPTt0(:,:,t), 'fro');
    nrm0(t) = trace(CGPTt0(:,:,t));
end
figure; plot(linspace(0, Ntime0*dt0, Ntime0), nrm0); title('Trace of time-dependent CGPT: theoretical value');

%%
% Reconstruct CGPT and show error
CGPTr = zeros(ord*2, ord*2, Ntime);
for t=1:Ntime
    MSR = data.MSR{t};
    % out{t} = P.reconstruct_CGPT(MSR, ord, 100000, 1e-10, symmode);
    toto = P.reconstruct_CGPT_analytic(MSR, ord);
    CGPTr(:,:,t) = toto.CGPT;
    %     % out{f}.res/norm(MSR,'fro')
    %     norm(M{f} - out{f}.CGPT, 'fro')
end

nrm1 = zeros(1,Ntime);
for t=1:Ntime
    nrm1(t) = trace(CGPTr(:,:,t));
    % nrm1(t) = norm(CGPTr(:,:,t), 'fro');
end
figure; plot(linspace(0, Ntime*dt, Ntime), nrm1); title('Trace of time-dependent CGPT: reconstruction');

cc
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

