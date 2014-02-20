function [X, F, Q] = target_model(dt, Ntime, X0, std_acc, std_acc_angl)
%[X, F, Q] = target_model(dt, Ntime, X0, std_acc, std_acc_angl)
% Calculate the matices used in system state equation of KF and simulate the target path
%
% Inputs:
% dt: time period
% Ntime: Ntime X dt = total time
% X0: initial state vector
% std_acc, std_acc_angl: standard devirations for acceleration vector and angular acceleration
%
% Outputs:
% X: simulated target path, a matrix of dimesion 5 X Ntime. X(1:2,:) is the velocity vector,
% X(3:4,:) the position vector, X(5,:) the angular position
% F: system state matrix 
% Q: system noise covariance matix

F = [1  0 0 0 0; 
     0  1 0 0 0; 
     dt 0 1 0 0; 
     0 dt 0 1 0; 
     0  0 0 0 1];

sqdt = dt^2/2;

Q0 = [dt 0  0    0    0; 
      0  dt 0    0    0; 
      0  0  sqdt 0    0; 
      0  0  0    sqdt 0; 
      0  0  0    0    dt];

% F = [1 0 0 0 0; 0 1 0 0 0; dt 0 1 0 0; 0 dt 0 1 0; 0 0 0 0 1];
% sqdt = dt^2/2;
% Q0 = [dt 0 0 0 0; 0 dt 0 0 0; 0 0 sqdt 0 0; 0 0 0 sqdt 0; 0 0 0 0 dt];

% We generate first the system state vectors
X = zeros(5, Ntime);

% for n=1:Ntime
%     acc = std_acc * randn(2,1);
%     acc_angl = std_acc_angl * randn;
% 
%     if n==1
%         X(:,n) = F*X0 + Q0 * [acc; acc; acc_angl];
%     else
%         X(:,n) = F*X(:,n-1) + Q0 * [acc; acc; acc_angl];
%     end
% end

X(:,1) = X0;

for n=2:Ntime
    acc = std_acc * randn(2,1);
    acc_angl = std_acc_angl * randn;

    X(:,n) = F*X(:,n-1) + Q0 * [acc; acc; acc_angl];
end

acc_I2 = std_acc^2 * eye(2);
Sigma = [acc_I2 acc_I2 [0; 0]; 
         acc_I2 acc_I2 [0; 0]; 
         0 0 0 0 std_acc_angl^2];
Q = Q0 * Sigma * Q0';
