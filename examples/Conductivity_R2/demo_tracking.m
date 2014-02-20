%% Demo of tracking
% This script shows the tracking of a mobile target using extended Kalman filter. 
% Reference: Tracking of a mobile target using Generalized Polarization Tensors, SIAM Journal on imaging science, 2013

%% Add path
clear all;
close all;
clc;
addpath('../../');

%% Definition of the small inclusion

%%
% Initialize an object of |C2boundary|

% B = shape.Ellipse(1,1/2,2^9);
% B = shape.Flower(1/2, 1/2, 2^9);
% B = shape.Triangle(1/2, pi*0.8, 2^10);
% B = shape.Rectangle(1, 1/2, 2^10);
% B = shape.Banana(5/2, 1, [0,10]', 0, 0, 2^10);
B = shape.Imgshape('~/Data/images/Letters/A.png', 2^10);

%%
% Make inclusion
D = B*2;
cnd = 3;
pmtt = 0;

%% Trajectory of target
% We generate the trajectory for data simulation

dt = .01; % time step
ntime = 2 / dt; 

std_acc = 2; % standard deviation of acceleration vector
std_acc_angl = 1; % standard deviation of angular acceleration

%%
% Initial state of the system vector X=[v, z, a] of dimension 5, with v, z the vector of
% velocity and position, and a the angle (rotation). The initial state is
% biased to test the robustess of Kalman Filter.

% initial state
X0 = [-1, 1, 0.5, -0.5, pi/2]'; % X0 = zeros(5,1);
[X, F, Q] = tracking.target_model_bm(dt, ntime, X0, std_acc, std_acc_angl);

%% Set up an environment for experience
% The sources/receptors are distributed on a circle whose center is closed
% to the mass of center of the inclusion, up to a small offset.

% Radius of the measurement circle
rad = 1.5*(max(sqrt(X(3,:).^2 + X(4,:).^2)) + D.diameter);

cfg = acq.Coincided([0,0]', rad, 20, [1, 2*pi, 2*pi], 0);
% cfg = acq.Coincided([0,0]', rad, 20, [5, 0.2*pi, 2*pi], 0);

%% Simulate the MSR data stream

freq = 0;
MSR = zeros(cfg.data_dim, ntime);

for t=1:ntime
    Dt = (D<X(5,t)) + X(3:4,t);
    P = PDE.Conductivity_R2(Dt, cnd, pmtt, cfg);        
    toto = P.data_simulation(freq);
    MSR(:, t) = reshape(toto.MSR{1}, [], 1);
end

%% Configure the extended Kalman Filter (EKF)
% define the function handle h and dh

ord = 2; % max order of CGPT

% CGPT of the object of reference
lambda = asymp.CGPT.lambda(cnd, pmtt, freq);
CGPT_std = asymp.CGPT.theoretical_CGPT(D, lambda, ord);

% Forward system matrices
Sm = P.make_matrix_A(cfg.all_src, cfg.center, ord); % sources
Rm = P.make_matrix_A(cfg.all_rcv, cfg.center, ord); % receivers

H = @(x)tracking.func_h(x, Sm, Rm, CGPT_std);
dH = @(x)tracking.func_dh(x, Sm, Rm, CGPT_std);

% Add white noise to data and compute the covariance matrix of noise
sigma = 0.2 * norm(MSR, 'fro')/sqrt(ntime * cfg.data_dim);
MSR_noisy = MSR + randn(size(MSR))*sigma; 
R = eye(cfg.data_dim) * sigma^2;

%% Apply EKF for target tracking

% guess of initial state vector: Exp[X_0]
Xinit = X0 + 0.1*randn(5,1); 
% Xinit = zeros(5,1);

% guess of initial covariance matrix: Cov[X_0] 
P0 = zeros(5); 

Xt = tracking.Extended_Kalman_Filter(MSR, ntime, F, Q, H, dH, 2*R, Xinit, P0);

%% show result

fig=figure; 
plot(cfg, 'LineWidth', 1); hold on; axis image;
plot((D<X(5, end))+X(3:4, end)); plot(X(3,end), X(4,end),'rx');
plot((D<X(5, 1))+X(3:4, 1)); plot(X(3,1), X(4,1),'ro');
traj0 = plot(X(3,:), X(4,:));

traj1 = plot(Xt(3,:), Xt(4,:),'r');
legend([traj0, traj1],'True trajectory','Estimation')

fig=figure();
vx0 = plot(X(3,:)); hold on
vx1 = plot(Xt(3,:), 'r');
legend([vx0, vx1],'True position in x axis','Estimation')

fig=figure();
vy0 = plot(X(4,:)); hold on
vy1 = plot(Xt(4,:), 'r');
legend([vy0, vy1],'True position in y axis','Estimation')

fig=figure();
phase0 = plot(X(5,:)); hold on
phase1 = plot(Xt(5,:), 'r');
legend([phase0, phase1],'True orientation','Estimation')

fig=figure();
vx0 = plot(X(1,:)); hold on
vx1 = plot(Xt(1,:), 'r');
legend([vx0, vx1],'True velocity in x axis','Estimation')

fig=figure();
vy0 = plot(X(2,:)); hold on
vy1 = plot(Xt(2,:), 'r');
legend([vy0, vy1],'True velocity in y axis','Estimation')

