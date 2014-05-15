%% Demo of the ElectricFish class
% This script shows how to use |PDE.ElectricFish| class for data simulation
% and reconstruction of GPTs.

%% Add path
clear all;
close all;
clc;
addpath('../../');

%% Definition of the small inclusion

delta = 1 ; % diameter of the standard shape
%%
% Make (multiple) inclusions

nbPoints = 2^10; % number of discritization points on the boundary

% Initialize an object of |shape.C2boundary|
B = shape.Flower(delta/2, delta/2, nbPoints, 5, 0.4, 0);
D{1}=(B<(0.3*pi))*0.2+0.02*[1,1]'; % rotation, scaling and translation to get the first inclusion

cnd = [10]; % conductivity values of two inclusions
pmtt = [0.1]; % permittivity values

% One can also add multiple inclusions
% B = shape.Flower(delta/2, delta/2, nbPoints, 5, 0.4, 0);
% D{1}=(B<(0.3*pi))*0.2+0.2*[1,1]'; % rotation, scaling and translation to get the first inclusion
% B = shape.Ellipse(delta,delta/2,nbPoints);
% D{2}=B*0.2+0.2*[-1,-1]'; % second inclusion
% 
% cnd = [10, 0.1]; % conductivity values of two inclusions
% pmtt = [0.1, 0.2]; % permittivity values

%% Set up an environment for experience
% The sources/receptors are distributed on a circle whose center is closed
% to the mass of center of the inclusion, up to a small offset.
                           
%% Definition of fish's body
% We define the body of the fish, which is also an object of the
% |C2boundary| class.
%%
% The fish swims on a measurement circle around a measurement center, which
% is close to the center of mass of the inclusion, up to a small offset.

mcenter = [0,0]'; % center of movement
mradius = D{1}.diameter*2.5; % radius of the measurement circle

%%
% Initialize the fish's body, which may have different shapes.

Omega = shape.Banana(mradius*delta*2.5, mradius*delta/10, [mradius, 0]', 0, 0, nbPoints/2); % Banana-shaped fish
                                                                                     
% Omega = shape.Ellipse(delta*2, delta/4, 2^9); % Elliptic fish 
% Omega = Omega<(1/2*pi);

%%
% Set the skin's impedence

impd = 0.001;

%% 
% The fish's receptors are distributed on the skin, and one can choose 
% activate receivers by giving their indexes.

idxRcv = 1:2:Omega.nbPoints; % This generates a equally distributed receptors
% idxRcv = (Omega.nbPoints/4):4:(3*Omega.nbPoints/4); 

%% Configuration of the acquisition and setting up an environment for experience
% By configuration of the acquisition we mean the positions of the fish's
% electric organ and receptors, which depend on the trajectory. This is
% set up by the |acq.Fish_circle| class, which
% is a circular trajectory. Remark that the radius of measurement circle
% has no effect for banana-shaped fish.

cfg = acq.Fish_circle(Omega, idxRcv, mcenter, mradius, 10, 2*pi, [], [], 0.5, impd);
% figure; plot(cfg);

%%
% To speed up the numerical simulationof the P1 boundary element method
% (BEM), we down sampling the fish's body with a factor:
stepBEM = 4; % down-sampling factor for the P1 basis

%%
% Initialize an environment by passing the fish, the inclusions, the
% configuration, and the physical constants etc.
P = PDE.ElectricFish(D, cnd, pmtt, cfg, stepBEM);
figure; plot(P, [], 'LineWidth', 1); axis image;

%% Simulation of data
% The output is a structure containing many fields, see the function
% data_simulation().

%%
% The fish uses a list of working frequencies to acquire data.
freqlist = linspace(100,200,1); % frequency values

tic
data = P.data_simulation(freqlist);
toc

%% 
% After simulation, we verify that the solutions are in $L^2_0$ space.
% Their boundary integrals should be zero. 
%% 
% Boundaray integral of P1 elements $\int_{\partial\Omega} \psi(x)ds(x)$
% and P0 elements $\int_{\partial\Omega} \phi(x)ds(x)$
fprintf('Verify that the solution of the foward system (functions phi and psi) have zero-mean:');
norm(Omega.sigma * data.fpsi{1}(:,1)) % fpsi is the value of the function $\psi(x)$
norm(D{1}.sigma * data.fphi{1}(:,1, 1)) % fphi is the value of the function $\phi(x)$

%% Plot the potential fields
% *Bugs*: some grid points are on the fish's skin which may cause numerical
% instability. 
%%
% Calculate the field and plot it. 
sidx = 1; % source index to be calculated
fidx = 1; % frequency index to be calculated
[F, F_bg, SX, SY] = P.calculate_field(fidx, sidx, [0,0]', 3, 100, data.fpsi_bg, data.fpsi, data.fphi);
P.plot_field(sidx, F, F_bg, SX, SY, 100, 0, '-g','LineWidth', 1.1);

%% Reconstruction of CGPT
% Finally from data acquired by fish, we reconstruct the CGPTs and compare
% them with the theoretical values. Since data are of multi-frequency, and
% the GPTs are frequency dependent, we choose the working frequency before
% the reconstruction

fidx = 1; % frequency index to be reconstructed
symmode = 0; % force the solution to be symmetric
ord = 2; % maximum order of the reconstruction

%%
% Compute first the theoretical value of CGPT
for f=1:length(freqlist)
    lambda = asymp.CGPT.lambda(cnd, pmtt, freqlist(f));
    CGPTD{f} = asymp.CGPT.theoretical_CGPT(D, lambda, ord);
end

%% 
% The simulated data was noiseless, we add white noise to data
nlvl = 0.1; % noise level
data = P.add_white_noise(data, nlvl);

%%
% Reconstruct the CGPT from (noisy) data.
MSR = data.MSR_noisy{fidx}; % MSR data matrix
Cur = data.Current_noisy{fidx}; % Surface current

% or noiseless data
% MSR = data.MSR{fidx}; 
% Cur = data.Current{fidx};

%% 
% Pass all these arguments to reconstruct_CGPT() and also the number of
% iteration, accuracy and the constraint of symmetry to lsqr solver. The
% reconstruction is a structure.
out = P.reconstruct_CGPT(MSR, Cur, ord, 10000, 1e-10, symmode); 

fprintf('\n************ Exp. 1 **************\n');
fprintf('Theoretical CGPT matrix at the frequency %f:\n', freqlist(fidx));
CGPT0 = CGPTD{fidx} % theoretical value

fprintf('Reconstructed CGPT matrix at the frequency %f:\n', freqlist(fidx));
CGPT1 = out.CGPT % reconstruction

fprintf('Relative error\n');
norm(CGPT1-CGPT0,'fro')/norm(CGPT0,'fro') % relative error

%% Reconstruction of PT from the post-processed SFR (PP_SFR) matrix
% We can also reconstruct directly the first order PT by solving a
% different linear system which relates the post-processed SFR matrix to PT
% matrix (2x2).

% PP_SFR = data.PP_SFR{fidx};
PP_SFR = data.PP_SFR_noisy{fidx}; % Post-processed SFR data matrix
Cur_bg = data.Current_bg; % Surface current

out = P.reconstruct_PT(PP_SFR, Cur_bg, 10000, 1e-10, symmode); 

fprintf('\n************ Exp. 2 **************\n');
fprintf('Theoretical PT matrix at the frequency %f:\n', freqlist(fidx));
PT0 = CGPTD{fidx}(1:2,1:2) % theoretical value (first order)

fprintf('Reconstructed PT matrix at the frequency %f:\n', freqlist(fidx));
PT1 = out.PT % reconstructed value

fprintf('Relative error\n');
norm(PT1-PT0,'fro')/norm(PT0,'fro')

%% Reconstruction of the imaginary part of PT
% In the previous reconstuction of PT, one needs to calculate the background field U to a high
% precision. This procedure is error prone and not always feasible in practice. By taking the
% imaginary part of the data, we remove the background field U and reconstruct only the imaginary
% part of the PT.

out = P.reconstruct_PT(imag(PP_SFR), Cur_bg, 10000, 1e-10, symmode);

fprintf('\n************ Exp. 3 **************\n');

fprintf('Theoretical PT matrix (imaginary part only) at the frequency %f:\n', freqlist(fidx));
PT0_imag = imag(PT0)

fprintf('Reconstructed PT matrix (without background field U, imaginary part only) at the frequency %f:\n', freqlist(fidx));
PT1_imag = real(out.PT)

fprintf('Relative error\n');
norm(PT1_imag - PT0_imag, 'fro') / norm(PT0_imag, 'fro')

%% Compare the theoretical dipolar expansion value against the post-processed SFR
% We compare the measured data with the theoretical computations using the
% dipolar expansion. The red curve corresponds to the measurements.
% The error can be important if the size of the target is big or if the
% center of measurement is not close to the target.

PT_target = {}; % PT of the inclusion

for f=1:length(freqlist)
    lambda = asymp.CGPT.lambda(cnd, pmtt, freqlist(f));
    PT_target{f} = asymp.CGPT.theoretical_CGPT(D, lambda, 1);

    DiExp = reshape(out.op.Lnsym(PT_target{f}, 'notransp'), cfg.Ns_total, cfg.Nr);
    PP_SFR = data.PP_SFR_noisy{f};

    figure; 
    subplot(121); plot(real(DiExp(1,:))); hold on; plot(real(PP_SFR(1,:)), 'r'); title('Real part')
    subplot(122); plot(imag(DiExp(1,:))); hold on; plot(imag(PP_SFR(1,:)), 'r'); title('Imaginary part');
end
