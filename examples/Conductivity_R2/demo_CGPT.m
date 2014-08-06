%% Demo of the Conductivity_R2 class
% This script shows how to compute the contracted GPT (CGPT) matrix of
% inclusion(s). It shows also the stability of the CGPT matrix with respect
% to the (global or local) perturbation of the shape.

%% Add path
clear all;
close all;
clc;
addpath('../../');

%% Definition of small inclusions

%%
% Initialize an object of |C2boundary|

B = shape.Ellipse(1,1/2,2^9);
% B = shape.Flower(1/2, 1/2, 2^10, 5, 0.4, 0);
% B = shape.Triangle(1/2, pi*0.8, 2^10);
% B = shape.Rectangle(1,1/2,2^10);
% Omega = shape.Banana(5/2,1,[0,10]',0,0,1024);
% B = shape.Imgshape('../images/Letters/R.png', 1, 1, 2^10);

%%
% Make multiple inclusions
D{1}=(B<(0.5*pi))*0.5+[1,1]'; % rotation, scaling and translation
D{2}=B*0.5+[-1,-1]';
figure; plot(D{1}); hold on; plot(D{2}); axis image;

%% Compute the CGPT matrix of multiple inclusions

%%
% Constants of conductivity, permittivity and contrast of each inclusion
cnd = [0.5, 5]; 
pmtt = [1, 1];
lambda = asymp.CGPT.lambda(cnd, pmtt, 0);
ord = 3; % maximum order

%%
% Numerical evaluation of CGPT matrix. There is no linear relationships
% between the GPTs of multiple objects and those of each single object

fprintf('CGPT matrix of multiple inclusions:\n');
M = asymp.CGPT.theoretical_CGPT(D, lambda, ord) % multiple inclusions


fprintf('CGPT matrix of the first inclusion:\n');
M1 = asymp.CGPT.theoretical_CGPT(D{1}, lambda(1), ord) % single inclusion D1

fprintf('CGPT matrix of the second inclusion:\n');
M2 = asymp.CGPT.theoretical_CGPT(D{2}, lambda(2), ord) % single inclusion D2

%%
% Compare the numerical evaluation and the theoretical formula on an
% ellipse

fprintf('CGPT matrix of an ellipse computed by numerical integration:\n');
M_numeric = asymp.CGPT.theoretical_CGPT(B, lambda(1), ord) % single inclusion by numerical evaluation

fprintf('CGPT matrix of an ellipse computed by theoretical formula:\n');
M_formula = asymp.CGPT.ellipsetensor(ord, 1, 1/2, cnd(1)) % single inclusion by formula
% M2 = asymp.CGPT.disktensor(ord, delta, cnd(1))

fprintf('Relative error:\n');
norm(M_numeric-M_formula,'fro')/norm(M_formula,'fro') % numerical error error


%% CGPT of perturbed shapes.
% CGPT is stable to the global and the local perturbations

lambda = asymp.CGPT.lambda(3);
B = shape.Imgshape('~/Data/images/Letters/R.png', 2^10);
figure; plot(B); axis image;

fprintf('CGPT matrix of a letter:\n');
M0 = asymp.CGPT.theoretical_CGPT(B, lambda, ord)

%%
% Global perturbation
Bg = B.global_perturbation(0.01, 10, 50);
figure; plot(Bg); axis image;

fprintf('CGPT matrix of the globally perturbed shape:\n');
Mg = asymp.CGPT.theoretical_CGPT(Bg, lambda, ord) 

fprintf('Relative error:\n');
norm(Mg-M0,'fro')/norm(M0, 'fro')

%%
% Local perturbation
Bl = B.local_perturbation(-0.3, 0.35, 0.02, 0.1*pi);
figure; plot(Bl); axis image;

fprintf('CGPT matrix of the locally perturbed shape:\n');
Ml = asymp.CGPT.theoretical_CGPT(Bl, lambda, ord) 

fprintf('Relative error:\n');
norm(Ml-M0,'fro')/norm(M0, 'fro')

%% 
% Damage by smoothing a segment of the curve

Bl = B.smooth(0.2, 0.5, 0.2); % We smooth the segment corresponding to the parameterization 2*pi*[0.5-0.2, 0.5+0.2] using a constant window of length 0.2*(length of curve)
figure; plot(Bl); axis image;

fprintf('CGPT matrix of the damaged shape:\n');
Ml = asymp.CGPT.theoretical_CGPT(Bl, lambda, ord) 

fprintf('Relative error:\n');
norm(Ml-M0,'fro')/norm(M0, 'fro')
