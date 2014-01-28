function [X, F, Q] = target_model_angle(dt, Ntime, X0, std_acc_angl)
%[X, F, Q] = target_model(dt, Ntime, X0, std_acc, std_acc_angl)
% Calculate the matices used in system state equation of KF and simulate the target path
% Inputs:
% dt: time period
% Ntime: Ntime X dt = total time
% X0: initial state vector
% std_acc, std_acc_angl: standard devirations for acceleration vector and angular acceleration
% Outputs:
% X: simulated target path, a matrix of dimesion 5 X Ntime. X(1:2,:) is the velocity vector,
% X(3:4,:) the position vector, X(5,:) the angular position
% F: system state matrix 
% Q: system noise covariance matix

F = 1;

% We generate first the system state vectors
X = zeros(1, Ntime);

for n=1:Ntime
    acc_angl = std_acc_angl * randn;

    if n==1
        X(:,n) = F*X0 + dt * acc_angl;
    else
        X(:,n) = F*X(:,n-1) + dt * acc_angl;
    end
end

Q = std_acc_angl^2 * dt^2;
