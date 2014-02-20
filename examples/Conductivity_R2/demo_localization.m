%% Demo of Localization
% This script shows a simple projection method for localization of
% inhomogenehous inclusions in Electro Impedance Tomography

%% Add path
clear all;
% close all;
clc;
addpath('../../');

%% Definition of small inclusions

%%
% Initialize an object of |C2boundary|

B = shape.Ellipse(1, 1/2, 2^9);
% B = shape.Flower(1/2, 1/2, 2^10);
% B = shape.Triangle(1/2, pi*0.8, 2^10);
% B = shape.Rectangle(1, 1/2, 2^10);
% B = shape.Banana(5/2, 1, [0,10]', 0, 0, 2^10);
% B = shape.Imgshape('../images/Letters/R.png', 2^10);

%%
% Make multiple inclusions
D{1} = B*0.2 + 3.5*[1,1]';
D{2} = B*0.2 + 3.5*[-1,-1]';
cnd = [10, 10]; 
pmtt = [5, 5];

%% Set up an environment for experience
% The sources/receptors are distributed on a circle whose center is closed
% to the mass of center of the inclusion, up to a small offset.

%%
% Make the acquisition configuration with the class |acq.Coincided|.
cfg = acq.Coincided([0,0]', 7, 100, [1, 2*pi, 2*pi], 0);
% cfg = acq.Coincided([0,0]', 3, 10, [5, 0.2*pi, 2*pi], 0);

P = PDE.Conductivity_R2(D, cnd, pmtt, cfg); 

fig=figure; plot(P, 'LineWidth', 1); axis image;

%% Simulation of the MSR data
% freqlist = linspace(0, 0.1, 1);
freqlist = 10;
tic
data = P.data_simulation(freqlist);
toc

%% Analysis of MSR matrix

%%
% SVD of the MSR matrix. We compute the real matrix in the complex case.
MSR = data.MSR{1};
[UU,SS,VV] = svd(real(MSR*MSR')); 

%%
% Determines how many significant non zeros singular values, which is also the number of
% inclusions. We suppose this is done correctly.
ss = sqrt(diag(SS));
figure; semilogy(ss);

N0 = 2*length(D); % Number of inclusions
U = UU(:, 1:N0);  % corresponded singular vectors

%% Localization by projection

Ns = 100; % number of searching points
b0 = 4; % searching box is [-b0, b0]^2
sx = linspace(-b0, b0, Ns);
sy = linspace(b0, -b0, Ns);
[Sx, Sy] = meshgrid(sx, sy);
Err = zeros(Ns);

for n=1:Ns
    for m=1:Ns
        % construct the test vector
        Z0 = [Sx(n,m); Sy(n,m)];        
        V = cfg.all_src - repmat(Z0, 1, cfg.Ns_total); V = V';
        nV = sqrt(V(:,1).^2+V(:,2).^2);
        
        % the vector must be normalized
        vec1 = V(:,1)./nV; vec2 = V(:,2)./nV;
        % vec1 = V(:,1); vec2 = V(:,2);
        
        % Projection onto the image space of singular vectors
        p1 = vec1 - U*pinv(U)*vec1;
        p2 = vec2 - U*pinv(U)*vec2;

        Err(n,m) = 1/(norm(p1)^2 + norm(p2)^2);
    end
end

%%
% Show the image
figure; imagesc(Err); axis image; colorbar();

