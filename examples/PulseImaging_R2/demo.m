%% Demo of the PulseImaging_R2 class
% This script shows how to use |PDE.PulseImaging| class for data simulation
% and reconstruction of time dependent CGPTs.

%% Add path
clear all;
close all;
addpath('../../');

%% Definition of the small inclusion

%%
% Make (multiple) inclusions

nbPoints = 2^10; % number of discritization points on the boundary

% Initialize an object of |shape.C2boundary|
B = shape.Ellipse(1,1/2,nbPoints);
% B = shape.Triangle(1/2, pi*0.8, nbPoints);
% B = shape.Rectangle(1, 1/2, nbPoints);
% B = shape.Flower(1/2,1/2,nbPoints); % flower
% B = shape.Imgshape('~/Data/images/Letters/A.png', nbPoints);

%%
% True shapes
% D{1}=(B<(0.3*pi))*0.2+0.2*[1,1]'; % rotation, scaling and translation to get the first inclusion
D{1}=(B<(0.3*pi))*0.5; % rotation, scaling and translation to get the first inclusion

cnd = [3]; % conductivity values of two inclusions
pmtt = [1]; % permittivity values

% One can also add multiple inclusions
% B = shape.Ellipse(delta,delta/2,nbPoints);
% D{2}=(B<(0.3*pi))*0.2+2*[-1,-1]'; % second inclusion

% cnd = [10, 5]; % conductivity values of two inclusions
% pmtt = [1, 2]; % permittivity values

%% Set up an environment for experience
% The sources/receptors are distributed on a circle whose center is closed
% to the mass of center of the inclusion, up to a small offset.

cfg = acq.Coincided([0,0]', 3, 10, [1, 2*pi, 2*pi], false, [1,-1]);
% cfg = acq.Coincided([0,0]', 3, 10, [1, 0.5*pi, 2*pi], false, [1,-1]);

%% Show the pulse waveform h and its Fourier transfom H
Ntime = 2^8; % number of time steps
Tmax0 = 8; % time duration of h at the scale 0
[waveform, dt, Tmax, freqform, df, Fmax] = tools.make_pulse(Tmax0, Ntime, 1);

figure; 
subplot(211); plot(linspace(0, Tmax, Ntime), waveform); title('h'); xlabel('time');

Nfreq = length(freqform);
subplot(212); plot(linspace(0, Fmax, Nfreq), abs(freqform)); title('H'); xlabel('frequency');

%%
% Initialize an environment by passing the fish, the inclusions, the
% configuration, and the physical constants etc.

P = PDE.PulseImaging_R2(D, cnd, pmtt, waveform, dt, cfg);
figure; plot(P, 'LineWidth', 1); axis image;

%% Simulation of data
% The output is a structure containing many fields, see the function
% data_simulation().
disp('Data simulation...');

tic
data = P.data_simulation(Ntime);
toc

%% Reconstruction of time-dependent CGPT
%%
% add white noise
nlvl = 0.1; 
data = P.add_white_noise(data, nlvl);

%%
% Compute theoretical values of time-dependent CGPTs
disp('Computation of theoretical time dependent CGPTs...');

ord = 4; % maximal order

tic
[CGPTt0, dt0, CGPTf0] = asymp.CGPT.theoretical_CGPT_time(D, cnd, pmtt, ord, freqform, df, Tmax, Ntime);
toc

%%
% Reconstruct CGPT
disp('Reconstruction of time dependent CGPTs from data...');

ord = 2; % maximum order of the reconstruction
symmode = 1; % Force the solution to be symmetric (for lsqr method only)

tic
out = P.reconstruct_CGPT_analytic(data.MSR, ord);
%out = P.reconstruct_CGPT(data.MSR, ord, 100000, 1e-10, symmode, 'pinv');
toc

CGPTtr = zeros(2*ord, 2*ord, Ntime);

for t=1:Ntime
    CGPTtr(:,:,t) = out.CGPT{t};
end

%%
% Comparaison and show error

rr = 1; cc = 1;
nrm0 = squeeze(CGPTt0(rr,cc,:));
nrm1 = squeeze(CGPTtr(rr,cc,:));

% nrm0 = zeros(1,Ntime);
% for t=1:Ntime
%      nrm0(t) = trace(CGPTt0(1:2*ord, 1:2*ord, t));
% end
% figure; plot(linspace(0, Ntime0*dt0, Ntime0), nrm0); 
% title('Trace of time-dependent CGPT: theoretical value');

% nrm1 = zeros(1,Ntime);
% for t=1:Ntime
%     nrm1(t) = trace(CGPTtr(1:2*ord, 1:2*ord, t));
%     % nrm1(t) = norm(CGPTtr(1:2*ord, 1:2*ord, t), 'fro');
% end
% figure; plot(linspace(0, Ntime*dt, Ntime), nrm1);
% title('Trace of time-dependent CGPT: reconstruction');

figure; 
plot(linspace(0, Tmax, Ntime), nrm1); hold on;  
plot(linspace(0, Tmax, Ntime), nrm0, 'r'); 
plot(linspace(0, Tmax, Ntime), 0.6237*waveform, 'g'); 
title('Trace of time-dependent CGPT: comparison');

disp('Relative error:')
norm(nrm0-nrm1)/norm(nrm0) % This error decays as 1/Ntime, 

xx
%% Plot the potential fields and make a movie
%%
% Calculate the field and plot it. 
sidx = 1; % source index to be calculated
[F, F_bg, SX, SY] = P.calculate_field(Ntime, sidx, [0,0]', 3, 200);

% P.plot_field(sidx, F{20}, F_bg{20}, SX, SY, 20); % , '-g','LineWidth', 1.5);

%%
% Make a movie

Img = {};
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

%% The cost of initializing these objects is about 0.6s
% tic
% for tt=1:100
%     cfg = acq.Coincided([0,0]', 3, 10, [1, 2*pi, 2*pi], false, [1,-1]);
%     
%     Ntime = 2^9; % time interval length
%     [waveform, dt, Tmax, freqform, df, Fmax] = tools.make_pulse(Ntime, 3, 1);
%         
%     P = PDE.PulseImaging_R2(B{1}, cnd, pmtt, waveform, dt, cfg);
% end
% toc
