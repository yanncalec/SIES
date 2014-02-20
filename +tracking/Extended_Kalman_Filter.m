function  X_tt = Extended_Kalman_Filter(Y, Ntime, F, Q, h, dh, R, X0, P0)
% X_tt = Extended_Kalman_Filter(Y, Ntime, F, Q, h, dh, R, X0, P0)
% Extended Kalman Filter
% Inputs:
% Y: data stream of dimension ? X Nt, with Nt >= Ntime. If a cell of
% matrices is given, it will be converted to a data stream matrix.
% Ntime: processing time
% F: system state matrix 
% Q: system noise covariance matix
% h, dh: observation function and its derivatives, function handles
% R: observation noise covariance matrix
% X0: guess for initial state
% P0: guess for covariance matrix of the initial state

X_tt = zeros(length(X0), Ntime);
%X_tt(:, 1) = X0;

% Q0 = zeros(size(Q,1), 1);
% R0 = zeros(size(R,1), 1);

if nargin <= 8 
    P0 = zeros(size(Q));
end

% if nargin <= 9
%     verbose = 0;
% end
% 
P_tt = P0;

for n=1:Ntime
    %     if verbose
    %         disp(['Tracking time ', num2str(n), ' of ',num2str(Ntime)]);
    %     end
    
    % Prediction
    if n==1
        X_tn = F * X0;
    else
        X_tn = F * X_tt(:, n-1);
    end
    
    yt = Y(:,n) - h(X_tn);
    H = dh(X_tn);
    
    P_tn = F * P_tt * F' + Q;
    % if n==1
    %     P_tn = F * P0 * F' + Q;
    % else
    %     P_tn = F * P_tt * F' + Q;
    % end

    % Update    
    St = H * P_tn * H' + R;
    Kt = P_tn * H' * inv(St);
    X_tt(:, n) = X_tn + Kt * yt;
    P_tt = P_tn - Kt * H * P_tn;
end
