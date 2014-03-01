%% Demo of dictionary matching with the Helmholtz class
% This script shows how to use |PDE.Helmhotz_R2| class for dictionary matching

%% Add path
clear all;
close all;
clc;
addpath('~/OOP/');

load ~/Data/smalldico.mat;
load ~/Data/MSR.mat;

ord = (size(Dico{1}.SCT{1},1)-1)/2;

%% Matching test
nlvl = 0.5; % level of noise

for idx=1:length(D)
    S = zeros(Nv, Nv, length(sfreq));
    
    for n=1:length(sfreq)
        Data{idx} = P.add_white_noise(Data{idx}, nlvl);
        MSR = Data{idx}.MSR_noisy{n};
        P.freq = sfreq(n);
        out = P.reconstruct_SCT_analytic(MSR, ord);
        S(:,:,n) = dico.Helmholtz.ShapeDescriptorSCT(out.SCT, Nv);
    end
    
    [err1, idx1, scl1] = dico.Helmholtz.SCT_matching(S, sfrange, Dico.SD_S, frange);

    fprintf('True shape: %s\n', D{idx}.name_str);
    fprintf('Identified shapes: %s, %s, %s\n', D{idx1(1)}.name_str, D{idx1(2)}.name_str, D{idx1(3)}.name_str);
    fprintf('Errors: %f, %f, %f\n', err1(idx1(1)), err1(idx1(2)), err1(idx1(3)));
    fprintf('Estimations of scaling: %f, %f, %f\n\n', scl1(idx1(1)), scl1(idx1(2)), scl1(idx1(3)));
end

%% Matching with MSR
nlvl = 0; % level of noise
n=1;
idx=1;
MSR = Data{idx}.MSR{n};
P.freq = sfreq(n);
out = P.reconstruct_SCT_analytic(MSR, ord);
S = dico.Helmholtz.ShapeDescriptorSCT(out.SCT, Nv);

figure; imagesc(flipud(S));
figure; imagesc(abs(MSR));

%% plot
rr = 1; cc = 1;
T1 = squeeze(S(rr,cc,:)); toto = Dico.SD_S{3}; T2=squeeze(toto(rr,cc,:));

figure; plot(T1);  figure; plot(T2); 

[serr, vscl, idx1] = dico.Helmholtz.scaling_lookup_table(T2, frange, T1, sfrange);
vscl(idx1)
figure; plot(serr)
% [err1, idx1, scl1] = dico.Helmholtz.SCT_matching(T1, sfrange, {T2}, frange)

fidx=109;
D1 = (D{idx}<rtn) * scl + trl;
tic
W1 = D1.SCT(ord, sfreq(fidx), pmtt_bg, pmeb_bg);
toc
S1 = dico.Helmholtz.ShapeDescriptorSCT(W1, Nv);
norm(S(:,:,fidx) - S1, 'fro')/norm(S1,'fro')

%%

sfrange = scanfreq_range;
frange = [freqlist(1), freqlist(end)];

scl = 1.5; trl = [-0.5,0.5]'; rtn = pi/3;
D1 = (D{3}<rtn) * scl + trl;
%figure; plot(D1);

D1.SCT(ord, )
dico.ShapeDescriptor();
load WB;
load WD;

%% Definition of the small inclusion

delta = 1 ; % diameter of the standard shape
%%
% Initialize an object of |C2boundary|

% B = shape.Ellipse(delta,delta/2,[0,0]',0,2^10); % 
B = shape.Flower(delta/2,delta/2,[0,0]',0,2^10,5,0.4,0.5); 
% B = shape.Triangle(delta/2, pi*0.8, 2^10);
% B = shape.Rectangle(delta,0.5*delta,[0,0]',0,2^10);
% % or load image from file
% B = shape.Imgshape('../images/Letters/R.png', delta, delta, 2^10);
% figure; plot(B); axis image;

%%
% Set the conductivity and the permittivity of the inclusion
B.cnd = 3; B.pmtt = 3; B.pmeb = 3;

%% Compute the SCT of a reference shape

%%
% % Working frequencies of sources and background
freqlist = linspace(0.1, 20, 256)*pi; % Frequency range: from 0.1*pi to 20*pi
pmtt_bg = 1; 
pmeb_bg = 1;

N0 = 50;
WB={}; SB={}; GB={};

parfor n=1:length(freqlist)
    fprintf('Processing the frequency %d...\n',n);    
    WB{n} = B.SCT(N0, freqlist(n), pmtt_bg, pmeb_bg); 
end

for n=1:length(freqlist)
    [SB{n}, GB{n}] = dico.ShapeDescriptorSCT(WB{n});
end

fname = 'WB.mat';
save(fname, '-v7.3');

%% SCT of the rescaled shape

% Scaling, translation and rotation parameters:
% scl = 1; trl = [-0.5,0.5]'; rtn = pi/3;
scl = 2; trl = [0,0]'; rtn = 0;

%%
% The true inclusion _D_ is the shape _B_ after some rigid transform and
% scaling. Apply these transforms:
D = (B < rtn) * scl + trl; % or D = B.t_s_r(trl, scl, rtn);
figure; plot(D); axis image;

freqlist1 = linspace(1, 5, 64)*pi;

N0 = 50;
WD={}; SD={}; GD={};

parfor n=1:length(freqlist1)
    fprintf('Processing the frequency %d...\n',n);    
    WD{n} = D.SCT(N0, freqlist1(n), pmtt_bg, pmeb_bg);
end

for n=1:length(freqlist1)
    [SD{n}, GD{n}] = dico.ShapeDescriptorSCT(WD{n});
end

fname = 'WD.mat';
save(fname, '-v7.3');

%% Matching of scaling by lookup table
BM=zeros(256,256,length(freqlist));
for n=1:length(freqlist)
    BM(:,:,n)=SB{n};
end
DM=zeros(256,256,length(freqlist1));
for n=1:length(freqlist1)
    DM(:,:,n)=SD{n};
end

rr = 128; cc = 128;
rr = 32; cc = 64;
rr = 128; cc = 32;

rr = 1; cc = 1;

BM1 = reshape(BM(rr, cc, :), 1, []);
DM1 = reshape(DM(rr, cc, :), 1, []);
% BM1 = BM(1:4:end, 1:4:end, :);
% DM1 = DM(1:4:end, 1:4:end, :);

BM1 = mean(mean(BM,1),2);
DM1 = mean(mean(DM,1),2);
BM1 = reshape(BM1, 1, []);
DM1 = reshape(DM1, 1, []);

frange=[freqlist(1), freqlist(end)];
sfrange=[freqlist1(1), freqlist1(end)];

[err, scl1] = scaling_lookup_table_1d(BM1, frange, DM1, sfrange);
scl1


UUB = B.SCT(N0, 20, pmtt_bg, pmeb_bg);
SSB = dico.ShapeDescriptorSCT(UUB);
UUD = D.SCT(N0, 10, pmtt_bg, pmeb_bg);
SSD = dico.ShapeDescriptorSCT(UUD);

figure; imagesc(SSB);
figure; imagesc(SSD);
norm(SSB-SSD,'fro')

UUD1 = D.SCT(N0, 11, pmtt_bg, pmeb_bg);
SSD1 = dico.ShapeDescriptorSCT(UUD1);
norm(SSB-SSD1,'fro')/norm(SSB,'fro')
figure; imagesc(SSB-SSD1); colorbar()
figure; imagesc(SSD1);

%% Scaling invariant

rr=32; cc=32;
fmin = 2*pi; fmax = 10*pi;
idx0=find((freqlist>=fmin).*(freqlist<=fmax));

%tt0 = smooth(reshape(BM(rr,cc,idx0),1,[]),3);
tt0 = reshape(BM(rr,cc,idx0),1,[]);
% figure; plot(tt./(1:length(tt)));

fidx = freqlist(idx0);
figure; plot(fidx,tt0); hold on;
[xmax, imax, xmin, imin] = tools.extrema(tt0);
plot(fidx(imax),xmax,'r*',fidx(imin),xmin,'g*')

%tt1 = smooth(reshape(DM(rr,cc,:),1,[]),3);
tt1 = reshape(DM(rr,cc,:),1,[]);
% figure; plot(tt./(1:length(tt)));
fidx1 = freqlist1*2;
figure; plot(fidx1, tt1); hold on;
[xmax, imax, xmin, imin] = tools.extrema(tt1);
plot(fidx1(imax), xmax,'r*', fidx1(imin), xmin,'g*')

figure; plot(freqlist(idx0(1:2:end))); hold on; plot(2*freqlist1, 'r')
%% cc;

idx2=idx0(1:2:end);
ttsm = smooth(tt0,3);
tt2=ttsm(1:2:end);
figure; plot(idx2, tt2); hold on;
[xmax, imax, xmin, imin] = tools.extrema(tt2);
plot(idx2(imax),xmax,'r*',idx2(imin),xmin,'g*')
figure; plot(tt2); hold on;
figure; plot(tt0); hold on;



tt=[];
for n=1:128
    %    SSF= flipud(SS(:,:,n));
    tt(n,:) = diag(flipud(SS(:,:,n)));
    %    tt(n,:) = diag(SS(:,:,n));
    % tt(n)=norm(SSF(:,:,n),'fro');
end
figure; imagesc(tt);
%tt = reshape(SS(rr,cc,1:128),1,[]);
figure; plot(diag(SS(:,:,50)));

for rr=1:256
    for cc=1:256
        toto(rr,cc) = min(SS(rr,cc,:));
    end
end

for rr=1:256
    for cc=1:256
        toto(rr,cc) = min(SS(rr,cc,:));
    end
end

figure; imagesc(toto); colorbar()
figure; imagesc(SB{1}); colorbar()

%% 
figure;
for n=n0:n1
    subplot(2,4,n); 
    imagesc(abs(WB{n})); colorbar();
    axis image;
end

figure;
for n=1:4
    subplot(2,2,n); 
    imagesc(abs(WB{n+9-1})); colorbar();
    axis image;
end

figure;
for n=1:4
    subplot(2,2,n); 
    imagesc(abs(SB{n+9-1})); colorbar();
    axis image;
end

figure;
for n=n0:n1
    subplot(2,4,n); 
    imagesc(abs(SB{n})); colorbar();
    axis image;
end

figure;
for n=n0:n1
    subplot(2,4,n); 
    imagesc(abs(GB{n})); colorbar();
    axis image;
end

% WD1 = D.SCT(10, pi, pmtt_bg, pmeb_bg); 
% figure; imagesc(abs(WD1)); colorbar();
% figure; imagesc(abs(WD1-WB1)); colorbar();
% norm(WD1-WB1)

%% xxx
N0 = 60;
% WB={}; SB={}; GB={};

n0=17; n1=20;

for n=n0:n1
    WB{n} = B.SCT(N0, n*pi, pmtt_bg, pmeb_bg); 
end

figure;
for n=1:4
    subplot(2,2,n); 
    imagesc(abs(WB{n+n0-1})); colorbar();
    axis image;
end

figure;
for n=1:4
    subplot(2,2,n); 
    imagesc(abs(SB{n+n0-1})); colorbar();
    axis image;
end


%% maximum order of the reconstruction
ord = floor((N0-1)/2); 

disp('Computing the SCT matrix of the reference shape...');
for n=1:length(freqlist)
    fprintf('Processing the frequency %d...\n',n);
    W0{n} = B.SCT(ord, freqlist(n), pmtt_bg, pmeb_bg); 
end

disp('Computing the SCT matrix of the measured shape...');
for n=1:length(freqlist)
    fprintf('Processing the frequency %d...\n',n);
    W1{n} = D.SCT(ord, freqlist(n), pmtt_bg, pmeb_bg); 
end

cc
%% cc

for n=1:length(freqlist)
    W0{n} = B.SCT(ord, freqlist(n), pmtt_bg, pmeb_bg); 
end
W1 = E.SCT(ord, freq, pmtt_bg, pmeb_bg); 
tic
toc

figure; imagesc(abs(W0)); colorbar()
figure; imagesc(abs(W1)); colorbar()

%% Shape descriptors
[S0, G0] = dico.ShapeDescriptorSCT(W0);
[S1, G1] = dico.ShapeDescriptorSCT(W1);
figure; imagesc(S0); colorbar()
figure; imagesc(S0-S1); colorbar()
% figure; imagesc(([S1 S1; S1 S1])); colorbar()

disp('Error of the shape descriptor:')
norm(S0-S1,'fro')/norm(S0,'fro')

% norm(G0-G1,'fro')/norm(G0,'fro')
% norm(W0-W1,'fro')/norm(W0,'fro')

% figure; imagesc((G0)); colorbar()
% figure; imagesc(G0-flipud(fliplr(G0))); colorbar()

%% Set up an environment for experience
% The sources/receptors are distributed on a circle whose center is closed
% to the mass of center of the inclusion, up to a small offset.

%%
% offset of the measurement center
vh = D.diameter/2*(rand(2,1)-0.5); 
% vh= [0,0]';

%%
% center of the measurement circle
mcenter = D.center_of_mass + vh; 
%%
% radius of the measurement circle
mradius = D.diameter*2; 

%%
% Make the acquisition configuration. With the class |acq.Planewave|, we make equally distributed
% sources/receivers on a circle of center |mcenter| with radius |mradius|.

N0 = 21; % Number of sources (and receivers)
cfg = acq.Planewave(N0, mcenter, mradius, N0, [1, 2*pi, 2*pi]);
% figure;plot(cfg); axis image;

%%

P = PDE.Helmholtz_R2(D, cfg, freq, pmtt_bg, pmeb_bg); 

figure; plot(P, 'LineWidth', 1); axis image;

% xs = cfg.src(1); Sx = linspace(-3,3,100); 
% F = PDE.Conductivity_R2.solve_forward(D, 0, xs, Sx, Sx);
% figure; imagesc(F); colorbar()

%% Simulation of the MSR data
% freqlist = linspace(0,1,10);
tic
data = P.data_simulation();
toc

%%
% % Calculate the field and plot it. 
% [F, F_bg, SX, SY, mask] = P.calculate_field(3, cfg.src(1), 2*D.diameter, 64);
% P.plot_field(F, F_bg, SX, SY, 10, 'LineWidth', 2);

%%
% add white noise
nlvl = 0.1; 
data = P.add_white_noise(data, nlvl);
MSR = data.MSR{1}; 
% MSR = data.MSR_noisy{1}; 

%% Reconstruct SCT and show error
out = P.reconstruct_SCT_analytic(MSR, ord);

% out1 = P.reconstruct_SCT(MSR, ord);
% norm(out0.SCT0-out1.SCT0,'fro')/norm(out1.SCT0)

% figure; imagesc(abs(out.SCT0)); colorbar() % Without post-processing
% figure; imagesc(abs(out.SCT)); colorbar() % With post-processing

disp('Reconstruction error without post-processing:')
norm(out.SCT0-W1,'fro')/norm(W1,'fro')
disp('Reconstruction error with post-processing:')
norm(out.SCT-W1,'fro')/norm(W1,'fro')


%%%%%%

%% load Data
load smalldico.mat;
load MSR.mat;

