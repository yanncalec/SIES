%% Demo of the Conductivity_R2 class
% This script shows how to use |PDE.Conductivity_R2| class for data simulation
% and reconstruction of GPTs. The reconstruction of first order GPT (PT) is
% very stable wrt both the noise level and the angle of view.

% The reconstruction of the first order GPT (or the PT) is VERY ROBUST to
% the noise and to the angle of view in the acquisition system . For
% example, with the following setting:
% 50 transmitters, radius of measurement circle = 10 X radius of the object,
% and accept only 10% of relative error of reconstruction, then
% with pi/32 of aperture angle, one can even go beyond 100% of noise. In fact, the
% performance depends on several parameters:
% 1. Usage of the constraint of symmetry in the reconstruction, which can
% greatly enhance the robustness
% 2. Regularity  of the boundary. On those defined by analytic expressions 
% (ellipse, flower, triangle etc, with the corners smoothed using function 
% C2boundary.rescale) the performance is very good. However, the performance deteriorates
% dramatically on some irregular shapes like some letters (A, for example),
% while on some other shapes (letter E, L, for example) it still works
% well. It should be also noticed that the boundary smoothed by the function 
% C2boundary.smooth have systematically very bad performance. There may 
% exists bugs in the related functions. 
%
% In the limited view setting, the reconstruction of high order (>=2) is
% EXTREMELY UNSTABLE. 
%

%% Add path
clear all;
close all;
% addpath('../../');

% Letters='ABCDEFGHIJKLMNOPQRSTUVWXYZ';

%% Definition of small inclusions

%%
% Initialize an object of |C2boundary|

% B = shape.Ellipse(1,1/2,2^10);
% B = shape.Flower(1/2, 1/2, 2^10);
B = shape.Triangle(1/2, pi/3, 2^10, 10);
% B = shape.Rectangle(1, 1/2, 2^10);
% B = shape.Banana(2, 1/4, [0,10]', 0, 0, 2^10); B = B-B.center_of_mass;
% B = shape.Imgshape('~/Data/images/Letters/A.png', 2^10, 10); 

% B = B.smooth(10);
% B = B.local_perturbation(-0.05, 0.3, 0.05);

figure; plot(B); axis image;

%%
% Make (multiple) inclusion(s)

D{1}=(B<(pi/3))*1.5 + 0.1*[1,1]';
% D{1} = B;
% D{2}=B*0.5 + 0.3*[-1,-1]';

cnd = 10*[1, 1]; 
pmtt = 1*[1, 1];

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
cfg = acq.Coincided([0,0]', 10, 50, [1, 1/16*pi, 2*pi], false, [1,-1], 0.01);  

% Single Dirac
% cfg = acq.Coincided([0,0]', 10, 50, [1, 1/32*pi, 2*pi], false); 

% Full view
% cfg = acq.Coincided([0,0]', 4, 50, [1, 2*pi, 2*pi], 0);

% Non equally distributed
% cfg = acq.Coincided([0,0]', 4, 10, [5, 0.2*pi, 2*pi], 0);

P = PDE.Conductivity_R2(D, cnd, pmtt, cfg); 

fig=figure; plot(P, 'LineWidth', 1); axis image;

%% Simulation of the MSR data

freqlist = linspace(0, 100*pi, 5); % List of working frequencies
%freqlist = 0;

data = P.data_simulation(freqlist);
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
nlvl = 0;
nbExp = 1;
out = {};

K = max(1, ord);

% Reconstruction
for n=1:nbExp
    data = P.add_white_noise(data, nlvl);
    out{n} = P.reconstruct_CGPT(data.MSR_noisy, K, 100000, 1e-10, symmode, 'lsqr');
    % out{n} = P.reconstruct_CGPT(data.MSR_noisy, K, 100000, 1e-10, symmode, 'pinv');
end

fprintf('Relative error between theoretical and reconstructed CGPT matrix at different frequencies:\n');

for f=1:length(freqlist)
    % out{f}.res/norm(MSR,'fro')
    err = zeros(nbExp,1);
    errsvd = zeros(nbExp,1);

    for n=1:nbExp
        toto = out{n}.CGPT{f}(1:2*ord, 1:2*ord);
        x0 = svd(M{f}); x1=svd(toto); 
        x0 = x0(1)/x0(2); x1 = x1(1)/x1(2);
        err(n) = (norm(M{f} - toto, 'fro'))/norm(M{f},'fro');
        errsvd(n) = (norm(x0-x1, 'fro'))/norm(x1,'fro');
    end
    
    fprintf('Frequency: %f, error: %f, error of sv %f\n', freqlist(f), mean(err), mean(errsvd));
    % toto
    % toto - out.CGPT{f}
end

%% Error of MSR
% fprintf('Relative error between MSR and MSR from reconstructed CGPT matrix at different frequencies:\n');

% for f=1:length(freqlist)
%     toto = 0;
%     for n=1:nbExp
%         toto = out{n}.res{f} + toto;
%     end
%     toto = toto / nbExp;
    
%     fprintf('Frequency: %f, error: %f\n', freqlist(f), toto/norm(data.MSR{f}, 'fro'));
%     % toto
%     % toto - out.CGPT{f}
% end

%% Invariants of PT

% fprintf('PT is invariant to translation:\n');
% for f=1:length(freqlist)
%     lambda = asymp.CGPT.lambda(cnd, pmtt, freqlist(f));
%     M1{f} = asymp.CGPT.theoretical_CGPT(B+2*[1,1]', lambda, ord);
%     M0{f} = asymp.CGPT.theoretical_CGPT(B, lambda, ord);
% 
%     err = norm(M0{f}-M1{f},'fro')/norm(M0{f},'fro');
%     fprintf('Relative error between original and transformed PT matrix:%f\n', err);
% 
% end
% 
