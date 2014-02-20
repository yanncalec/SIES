function [X, F, Q] = target_model_bm(dt, Ntime, X0, std_acc, std_acc_angl)
% [X, F, Q] = target_model_bm(dt, Ntime, X0, std_acc, std_acc_angl)
% Calculate the matrices used in system state equation of KF and simulate
% the target path. The model is based on the brownian motion.
% 
% Inputs:
% dt: time step
% Ntime: number of time steps, Ntime X dt = total time
% X0: initial system state vector
% std_acc, std_acc_angl: standard devirations (scalar) for acceleration vector and angular acceleration
%
% Outputs:
% X: target path of simulation, a matrix of dimension 5 X Ntime. X(1:2,:) is the velocity vector,
% X(3:4,:) the position vector, X(5,:) the angular position
% F: system state matrix 
% Q: system noise covariance matix
%
% Reference [1]: Tracking of a mobile target using GPT


% Matrix of system state equation. Eq (25) of Ref[1]
F = [1  0 0 0 0; 
     0  1 0 0 0; 
     dt 0 1 0 0; 
     0 dt 0 1 0; 
     0  0 0 0 1];

% new model
b1 = std_acc^2 * eye(2);
b2 = std_acc^2 * dt * eye(2) /2;
b3 = std_acc^2 * dt^2 * eye(2) /3;
z1 = [0 0];

Q = dt*[b1 b2 z1';
    b2 b3 z1';
    z1 z1 std_acc_angl^2];

% [U,S,V] = svd(Q);
% Q0 = U*sqrt(S);
Q0 = sqrtm(Q);

% We generate first the system state vectors
X = zeros(5, Ntime);

X(:,1) = X0;
for n=2:Ntime
    U = Q0 * randn(5,1);    
    X(:,n) = F*X(:,n-1) + U;
end

% for n=1:Ntime
%     U = Q0 * randn(5,1);
%     
%     if n==1
%         X(:,n) = F*X0 + U;
%     else
%         X(:,n) = F*X(:,n-1) + U;
%     end
% end
