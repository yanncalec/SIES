%% Demo of the Conductivity_R2 class
% This script shows how to use |PDE.Conductivity_R2| class for data simulation
% and reconstruction of GPTs. The reconstruction of first order GPT (PT) is
% very stable wrt both the noise level and the angle of view.

% The reconstruction of the first order GPT (or the PT) is VERY ROBUST to
% the noise and to the angle of view in the acquisition system. For
% example, with the following setting:
% 50 transmitters, radius of measurement circle = 4 X radius of the object,
% and accept only 10% of relative error of reconstruction, then
% with 0.125pi of aperture angle, one can go up to 10% of noise
% with 0.25pi  of aperture angle, one can go up to 50% of noise
% with 0.5pi   of aperture angle, one can go up to 150% of noise
%
% In the limited view setting, the reconstruction of high order (>=2) is
% EXTREMELY UNSTABLE. 

%% Add path
clear all;
% close all;
% clc;
addpath('../../');
%% Definition of small inclusions

%%
% Initialize an object of |C2boundary|

% B = shape.Ellipse(1,1/2,2^9);
B = shape.Flower(1/2, 1/2, 2^9);
% B = shape.Triangle(1/2, pi*0.8, 2^10);
% B = shape.Rectangle(1, 1/2, 2^10);
% B = shape.Banana(5/2, 1, [0,10]', 0, 0, 2^10);
% B = shape.Imgshape('~/Data/images/Letters/R.png', 2^10);

fig=figure; 
plot(B); axis image;
saveas(fig,'~/Flower.eps','psc2');
%%
% Make (multiple) inclusion(s)
D{1} = B;
% D{1}=(B<(0.2*pi))*0.5 + 0.25*[1,1]';
% D{2}=B*0.5 + 0.3*[-1,-1]';
cnd = 5*[1, 1]; 
pmtt = 2*[1, 1];

% D{1}=B;
% % D{1}=(B<(0.5*pi))*0.5+[1,1]';
% cnd = [0.5]; 
% pmtt = [0];

%% Set up an environment for experience
% The sources/receptors are distributed on a circle whose center is closed
% to the mass of center of the inclusion, up to a small offset.

%%
% Make the acquisition configuration with the class |acq.Coincided|.

% limited angle of view

% Neutrality: surprisingly, this has a better conditionning
cfg = acq.Coincided([0,0]', 4, 50, [1, 1/8*pi, 2*pi], false, [1,-1], 0.01);  

% Single Dirac
%cfg = acq.Coincided([0,0]', 4, 50, [1, 0.25*pi, 2*pi], false); 

% Full view
% cfg = acq.Coincided([0,0]', 4, 50, [1, 2*pi, 2*pi], 0);

% Non equally distributed
% cfg = acq.Coincided([0,0]', 4, 10, [5, 0.2*pi, 2*pi], 0);

P = PDE.Conductivity_R2(D, cnd, pmtt, cfg); 

fig=figure; plot(P, 'LineWidth', 1); axis image;

%% Simulation of the MSR data

freqlist = linspace(0, 2*pi, 1); % List of working frequencies

tic
data = P.data_simulation(freqlist);
toc

%%
% Calculate and plot potiential fields

% sidx = 1;
% [F, F_bg, SX, SY, mask] = P.calculate_field(0.01, sidx, [0,0]', 6, 100);
% P.plot_field(sidx, F, F_bg, SX, SY, 100);

%% Reconstruction of CGPT
%%
% maximum order of the reconstruction
ord = 1;
symmode = 1;

%%
% Compute first the theoretical value of CGPT
for f=1:length(freqlist)
    lambda = asymp.CGPT.lambda(cnd, pmtt, freqlist(f));
    M{f} = asymp.CGPT.theoretical_CGPT(D, lambda, ord);
end

% % Verify that the matrix of CGPTs M is symmetric:
% sidx = 1
% fprintf('M-M^t at the frequency %f:\n', freqlist(sidx));
% M{sidx}-M{sidx}.'

%%
% Reconstruct CGPT and show error
% out = P.reconstruct_CGPT_analytic(data.MSR_noisy, ord);

% add white noise
nlvl = 0.2;
nbExp = 100;
out = {};

K = max(1, ord);

% Reconstruction
for n=1:nbExp
    data = P.add_white_noise(data, nlvl);
    out{n} = P.reconstruct_CGPT(data.MSR_noisy, K, 100000, 1e-10, symmode, 'pinv');
end

fprintf('Relative error between theoretical and reconstructed CGPT matrix at different frequencies:\n');

for f=1:length(freqlist)
    % out{f}.res/norm(MSR,'fro')
    err = zeros(nbExp,1);

    for n=1:nbExp
        toto = out{n}.CGPT{f}(1:2*ord, 1:2*ord);
        err(n) = (norm(M{f} - toto, 'fro'))/norm(M{f},'fro');
    end
    
    fprintf('Frequency: %f, error: %f\n', freqlist(f), mean(err));
    % toto
    % toto - out.CGPT{f}
end

fprintf('Relative error between MSR and MSR from reconstructed CGPT matrix at different frequencies:\n');

for f=1:length(freqlist)
    toto = 0;
    for n=1:nbExp
        toto = out{n}.res{f} + toto;
    end
    toto = toto / nbExp;
    
    fprintf('Frequency: %f, error: %f\n', freqlist(f), toto/norm(data.MSR{f}, 'fro'));
    % toto
    % toto - out.CGPT{f}
end


