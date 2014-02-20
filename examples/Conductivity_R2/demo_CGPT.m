%% Demo of the Conductivity_R2 class
% This script shows how to compute the contracted GPT (CGPT) matrix of multiple inclusions.

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

M = asymp.CGPT.theoretical_CGPT(D, lambda, ord) % multiple inclusions
M1 = asymp.CGPT.theoretical_CGPT(D{1}, lambda(1), ord) % single inclusion D1
M2 = asymp.CGPT.theoretical_CGPT(D{2}, lambda(2), ord) % single inclusion D2

%%
% Compare the numerical evaluation and the theoretical formula on an
% ellipse

M1 = asymp.CGPT.theoretical_CGPT(B, lambda(1), ord) % single inclusion by numerical evaluation
M2 = asymp.CGPT.ellipsetensor(ord, 1, 1/2, cnd(1)) % single inclusion by formula
% M2 = asymp.CGPT.disktensor(ord, delta, cnd(1))

norm(M1-M2,'fro')/norm(M2,'fro') % numerical error error


