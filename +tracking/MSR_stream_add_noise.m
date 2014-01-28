function [MSR, R] = MSR_stream_add_noise(MSR0, Ns, Nr, epsilon, symm)
% [MSR, R] = MSR_stream_add_noise(MSR0, Ns, Nr, epsilon)
% Add noise to a MSR stream
% Inputs:
% MSR0: noise free MSR stream of dimension (Ns X Nr) X Ntime
% Ns, Nr: dimension of a single MSR matrix (one frame in the time sequence)
% epsilon: percentage of noise
% Outputs:
% MSR: noisy MSR data
% R: the covariance matrix of the additive white noise

if nargin<5
	symm=false;
end
Ntime = size(MSR0, 2);
v = zeros(Ntime,1);

for n=1:Ntime
    toto = MSR0(:,n);
    v(n) = max(toto(:)) - min(toto(:));
end

vm  = mean(v); %mean value of MSR range

% Add noise
MSR = MSR0 + randn(size(MSR0)) * epsilon * vm;

if symm
	for n=1:Ntime
		toto = reshape(MSR(:, n), Ns, Nr);
		MSR(:, n) = reshape((toto + toto')/2, [], 1); % make symmetry
	end
end

% Covariance matrix of additive white noise
R = eye(Ns*Nr) * (epsilon*vm)^2;
