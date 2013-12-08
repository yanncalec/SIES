%% Demo of the Conductivity_R2 class
% This script shows how to compute the contracted GPT (CGPT) matrix of multiple inclusions.

%% Add path
clear all;
close all;
clc;
addpath('../../');

%% Definition of the small inclusion

delta = 1 ; % diameter of the standard shape
%%
% Initialize an object of |C2boundary|

B = shape.Ellipse(delta,delta/2,2^10);
% B = shape.Flower(delta/2, delta/2, 2^10, 5, 0.4, 0);
% B = shape.Triangle(delta/2, pi*0.8, 2^10);
% B = shape.Rectangle(delta,0.5*delta,2^10);

% % or load image from file
% B = shape.Imgshape('../images/Letters/R.png', delta, delta, 2^10);

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
ord = 5; % maximum order

%%
% Numerical evaluation of CGPT matrix
M = asymp.CGPT.theoretical_CGPT(D, lambda, ord); % multiple inclusions
M1 = asymp.CGPT.theoretical_CGPT(B, lambda(1), ord); % single inclusion
M2 = asymp.CGPT.ellipsetensor(ord, delta, delta/2, cnd(1)); % single inclusion by formula
% M2 = asymp.CGPT.disktensor(ord, delta, cnd(1))

norm(M1-M2,'fro')/norm(M2,'fro') % numerical error error


