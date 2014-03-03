%% Demo of the wavelet.OrthoWvl class
% This script shows how to use |wavelet.OrhoWvl| class for representation of features and how the
% imaging algorithm works.

%% Add path
clear all;
close all;
clc;
addpath('../../');

%% Definition of the small inclusion

%%
% Initialize an object of |C2boundary|

D = shape.Flower(1/2, 1/2, 2^10);
% D = shape.Ellipse(1,1,2^9);

cnd = 3; pmtt = 0;
%% Construction of Daubechies wavelet basis

%% 
% Width of the square ROI
ROI_width = D.diameter*1.2;

%%
% Resolution of the ROI
Rsl = ROI_width/32;

%% 
% Construction of the wavelet object

WF=wavelet.OrthoWvl('Daubechies', 6, D.center_of_mass, ROI_width, Rsl);

%%
% Number of scales for detail space. If 0, then only approximation space is used.
nbScl = 1;

%%
% Compute the wavelet coefficients (WPT)

disp('Computing theoretical WPT...');
tic
WPT = WF.WPT_std(D, cnd, nbScl);
toc

%% Plot the coefficient matrix (In the logarithmic scale)

X = WPT.AA; % The AA (approximation-approximation) coefficients

% X = WPT.DD{3,3}{1,1}; % The DD (detail-detail) coefficients

figure; imagesc(log10(abs(X))); axis image; colorbar();

%%
% Sorted AA coefficients (absolute value)
figure; semilogy(sort(abs(X(:)),'descend')); axis tight;
xlim([-0.1, 3]*10^5);

%%
% After keeping only the largest 1% of coefficients
Nz = ceil(0.01*numel(X));

[X1, Err1,~] = tools.Best_Napprx(X, Nz); 

Err1/norm(X,'fro') % Nonlinear approximation error

figure; imagesc(log10(abs(X1))); axis image; colorbar();

%% Localization property of the coefficients.
% We show the interaction of one wavelet (one row in the coefficient matrix) intersecting the
% boundary of the inclusion with all other wavelets. The interaction is localized, means that only
% the wavlets of the neighborhood have significant contribution. This property is commun for the AA
% and DD coefficients. The localization property becomes less significant when the conductivity
% constant of the inclusion increases (eg cnd = 50).

dim = WF.ApprxSpace{1}.dim; % Dimension of the coefficient matrix

% index of the wavelet
idx = 485; % For a flower
% idx = 347; % For a disk
Y0 = full(X(idx, :));

figure; plot(abs(Y0)); xlim([idx-150, idx+150]);
figure; imagesc(reshape(Y0, dim)); axis image; colorbar()
figure; imagesc(reshape(abs(Y0), dim)); axis image; colorbar()

%% Imaging from the wavelet coefficients

%%
% Imaging by the diagonal
Xim = WPT.diag_imaging(X);
figure; imagesc(abs(Xim)); axis tight; colorbar; title('Imaging of X');

%%
% Imaging by the maximum interaction. This imaging algorithm has better resolution than the
% previous one.
Xim = WPT.max_imaging(X);
figure; imagesc(abs(Xim)); axis tight; colorbar; title('Imaging of X');
% figure; imagesc(log10(abs(Xim))); axis tight; colorbar; title('Imaging of X in log');
