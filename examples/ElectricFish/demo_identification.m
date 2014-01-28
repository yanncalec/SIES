%% Demo of the shape identification using ElectricFish class
% This script shows the whole procedure of shape identification in a dictionary using |PDE.ElectricFish| class

%% Add path
clear all;
close all;
clc;
addpath('../../');

%% Definition of the small inclusion

delta = 1 ; % diameter of the standard shape
%%
% Initialize an object of |C2boundary|
nbPoints = 2^9;
Dico{1} = shape.Flower(delta/2, delta/2, nbPoints, 5, 0.4, 0.9);
Dico{2} = shape.Ellipse(delta,delta/2,nbPoints);
Dico{3} = shape.Triangle(delta/2, pi*0.8, nbPoints);
Dico{4} = shape.Rectangle(delta,0.5*delta,nbPoints);
% Dico{5} = shape.Imgshape('~/Data/images/Letters/R.png', delta, delta, nbPoints);

%%
% Plot the dictionary
figure;
subplot(221); plot(Dico{1});
subplot(222); plot(Dico{2});
subplot(223); plot(Dico{3});
subplot(224); plot(Dico{4});
% axis image;

%%
% Make an inclusion

D{1}=(Dico{1}<(0.3*pi))*0.5+[.1,-.1]';
cnd = [10]; 
pmtt = [0.1];

%% Set up an environment for experience
% The sources/receptors are distributed on a circle whose center is closed
% to the mass of center of the inclusion, up to a small offset.
                           
%% Definition of fish's body
% We define the body of the fish, which is also an object of the
% |C2boundary| class.
%%
% The fish swims on a measurement circle around a measurement center, which
% is close to the center of mass of the inclusion, up to a small offset.

mcenter = [0,0]';
mradius = D{1}.diameter*1.5; % radius of the measurement circle

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
idxRcv = 1:1:Omega.nbPoints; % This generates a equally distributed receptors
% idxRcv = (Omega.nbPoints/4):4:(3*Omega.nbPoints/4);

%% Configuration of the acquisition and setting up an environment for experience
% By configuration of the acquisition we mean the positions of the fish's
% electric organ and receptors, which depends on the trajectory. This is
% represented for example by an object of the |acq.Fish_circle| class, which
% is a circular trajectory. Remark that the radius of measurement circle
% has no effect for banana-shaped fish.
cfg = acq.Fish_circle(Omega, idxRcv, mcenter, mradius, 10, 2*pi, [], [], 0.5, impd);
% figure; plot(cfg);

%%
% To speed up the numerical simulationof the P1 boundary element method
% (BEM), we down sampling the fish's body with a factor
stepBEM = 4; % down-sampling factor for the P1 basis

%%
% Initialize an environment by passing the fish, the inclusion, the
% configuration, and the working frequencies etc.
P = PDE.ElectricFish(D, cnd, pmtt, cfg, stepBEM);
figure; plot(P, [], 'LineWidth', 1); axis image;

%% Simulation of data
% The output is a structure containing many fields, see the function
% data_simulation().

%%
% The fish uses a list of working frequencies to acquire data.
% freqlist = [0.01]; % one frequency            
freqlist = linspace(100,200,20); % multi-frequency

tic
data = P.data_simulation(freqlist);
toc

%% Construct the dictionary of shape descriptor and PT

ord = 3; % order of the dictionary

% Multi-frequency dictionary of CGPT, PT and imaginary part of PT
CGPT_Dico={}; 
PT_Dico = {};
PT_Dico_imag = {};

SD1_Dico={}; SD2_Dico={}; % Multi-frequency dictionary of the shape-descriptor

for f=1:length(freqlist)
    lambda = asymp.CGPT.lambda(cnd, pmtt, freqlist(f));
    toto = {};
    toto1 = {}; toto2={};

    for n=1:length(Dico)
        % CGPT ls computed separately for each shape
        toto{n} = asymp.CGPT.theoretical_CGPT(Dico{n}, lambda, ord);
        [U1{n}, U2{n}] = dico.CGPT.ShapeDescriptor_CGPT(toto{n});
        toto1{n} = toto{n}(1:2,1:2);
        toto2{n} = imag(toto1{n});
    end
    
    CGPT_Dico{f} = toto;
    SD1_Dico{f} = U1; SD2_Dico{f} = U2;

    PT_Dico{f} = toto1;
    PT_Dico_imag{f} = toto2;
end

%% Reconstruction of CGPT and shape identification
%% 
% Add white noise to data
nlvl = 0.5; % noise level
data = P.add_white_noise(data, nlvl);

%%
% Reconstruct the CGPT of a fixed frequency
fidx = 5; % frequency index for reconstruction of CGPT
symmode = 1; % force the solution to be symmetric

% MSR = data.MSR{fidx}; 
MSR = data.MSR_noisy{fidx}; % MSR data matrix
                            
%Cur = data.Current_noisy{fidx};
Cur = data.Current{fidx};

out = P.reconstruct_CGPT(MSR, Cur, ord, 10000, 1e-10, symmode); 

%%
% Compare the reconstructed shape descriptor with the dictionary

[I1, I2] = dico.CGPT.ShapeDescriptor_CGPT(out.CGPT); % make shape descriptor
cord = 2; % order of SD for matching

[err, idx] = dico.CGPT.SD_Matching(I1, I2, SD1_Dico{fidx}, SD2_Dico{fidx}, cord);

for n=1:cord
    fprintf('Identified shape by the shape descriptor of order %d is: %s\n', n, ...
            Dico{idx{n}(1)}.name_str);
end

%% Shape identification using PT 
%%
% Compute the multifrequency dictionary of PT

CGPT = {}; % Multi-frequency CGPT of reconstruction
PT = {}; % PT of reconstruction
PT_imag = {}; % imaginary PT of reconstruction

for f=1:length(freqlist)
    MSR = data.MSR_noisy{f}; % MSR data matrix    
                             
    Cur = data.Current_noisy{f};    
    % Cur = data.Current{f};    

    out = P.reconstruct_CGPT(MSR, Cur, 1, 10000, 1e-10, symmode); % reconstruct only the first order
    CGPT{f} = out.CGPT;

    % PP_SFR = data.PP_SFR{f};
    PP_SFR = data.PP_SFR_noisy{f};
    % Cur_bg = data.Current_bg;
    Cur_bg = data.Current_bg_noisy;
    
    out = P.reconstruct_PT(PP_SFR, Cur_bg, 10000, 1e-10, symmode); 
    PT{f} = out.PT;

    out = P.reconstruct_PT(imag(PP_SFR), Cur_bg, 10000, 1e-10, symmode); 
    PT_imag{f} = real(out.PT);
end

% Multi-frequency matching
[err, idx] = dico.CGPT.MF_PT_Matching(CGPT, PT_Dico);
fprintf('Identified shape by the multifrequency CGPT is: %s\n', Dico{idx(1)}.name_str);

[err, idx, PT_inv, PT_Dico_inv] = dico.CGPT.MF_PT_Matching(PT, PT_Dico);
fprintf('Identified shape by the multifrequency PT is: %s\n', Dico{idx(1)}.name_str);

[err, idx] = dico.CGPT.MF_PT_Matching(PT_imag, PT_Dico_imag);
fprintf('Identified shape by the multifrequency imag PT is: %s\n', Dico{idx(1)}.name_str);

