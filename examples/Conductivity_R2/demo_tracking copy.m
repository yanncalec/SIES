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
ntime = 1 / dt; 

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
data = tracking.MSR_simulation_stream(cfg, D, X(3:4,:), X(5,:));
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

% Use a renormalization constant to compensate numerical problem
% nrmlcst = 1/(max(MSR(:))-min(MSR(:)));
nrmlcst = 1; %1/mean(abs(MSR(:)));

H = @(x)tracking.func_h(x, Sm, Rm, CGPT_std, nrmlcst);
dH = @(x)tracking.func_dh(x, Sm, Rm, CGPT_std, nrmlcst);

% Add white noise to data and compute the covariance matrix of noise

% sigma = 0.1 * norm(MSR, 'fro')/sqrt(ntime * cfg.data_dim);
% MSR_noisy = MSR + randn(size(MSR))*sigma; 
% R = eye(cfg.data_dim) * sigma^2;

MSR = MSR0 * nrmlcst;

[MSR_noisy, R] = tracking.data_stream_add_noise(MSR, cfg.Ns_total, cfg.Nr_total, 0.1);

% data = P.add_white_noise(data, 0.1);
% R = eye(cfg.data_dim) * mean(data.sigma)^2*100;
% MSR_noisy = data.MSR_noisy;

%% Apply EKF for target tracking

% guess of initial state vector: Exp[X_0]
% Xinit = [0.5,0.5,1.0,1.0,0]';
% Xinit = X0 + 0*randn(5,1); 
Xinit = zeros(5,1);

% guess of initial covariance matrix: Cov[X_0] 
P0 = zeros(5); 
% P0 = eye(5);

Xt = tracking.Extended_Kalman_Filter(MSR, ntime, F, Q, H, dH, 2*R, Xinit, P0);

% show result
close all;
fpath='figures/lim_aov_1';

fig=figure; 
plot(P, 'LineWidth', 1); hold on; axis image;
traj0 = plot(X(3,:), X(4,:));
plot(D+X(3:4, end)); plot(X(3,end), X(4,end),'rx');

traj1 = plot(Xt(3,:), Xt(4,:),'r');
legend([traj0, traj1],'True trajectory','Estimation')

fig=figure();
vx0 = plot(X(3,:)); hold on
vx1 = plot(Xt(3,:), 'r');
legend([vx0, vx1],'True position in x axis','Estimation')
% saveas(fig,[fpath,'/pos_x.eps'], 'psc2');

fig=figure();
vy0 = plot(X(4,:)); hold on
vy1 = plot(Xt(4,:), 'r');
legend([vy0, vy1],'True position in y axis','Estimation')
% saveas(fig,[fpath,'/pos_y.eps'], 'psc2');

fig=figure();
phase0 = plot(X(5,:)); hold on
phase1 = plot(Xt(5,:), 'r');
legend([phase0, phase1],'True orientation','Estimation')
% saveas(fig,[fpath,'/orientation.eps'], 'psc2');

fig=figure();
vx0 = plot(X(1,:)); hold on
vx1 = plot(Xt(1,:), 'r');
legend([vx0, vx1],'True velocity in x axis','Estimation')
% saveas(fig,[fpath,'/v_x.eps'], 'psc2');

fig=figure();
vy0 = plot(X(2,:)); hold on
vy1 = plot(Xt(2,:), 'r');
legend([vy0, vy1],'True velocity in y axis','Estimation')
% saveas(fig,[fpath,'/v_y.eps'], 'psc2');


% plot(Xs(1,:), Xs(2,:), 'bx'); hold on; axis image;
% plot(X(3,1), X(4,1),'ro'); % initial position
% traj1 = plot(Xt(3,:), Xt(4,:),'r');
% legend([traj0, traj1],'True trajectory','Estimation')


% fig=figure();
% plot(Xs(1,:), Xs(2,:), 'bx'); hold on; axis image;
% plot(X(3,1), X(4,1),'ro');
% traj0 = plot(X(3,:), X(4,:));
% 
% plot(D_L{1}(1,:), D_L{1}(2,:), 'g');
% n=ceil(ntime/2);
% plot(X(3,n), X(4,n),'rx');
% plot(D_L{n}(1,:), D_L{n}(2,:), 'g');
% plot(X(3,end), X(4,end),'rx');
% plot(D_L{end}(1,:), D_L{end}(2,:), 'g');
% %saveas(fig,[fpath,'/tracking.eps'], 'psc2');
% 

