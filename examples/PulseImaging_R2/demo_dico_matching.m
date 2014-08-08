%% Make a dictionary of shape descriptors

%%
clear all;
close all;
addpath('~/SIES');

% Path to image file of letters
imagepath = '~/Data/images/Letters';

%% Parameters for the shape
nbPoints = 2^10; % Number of boundary points for discretization
delta = 1; % standard size

%%
% All shapes have the same conductivity and permittivity values
cnd = 5;
pmtt = 2;

%% Make dictionary

disp('Construction of the dictionary...');

B{1} = shape.Ellipse(delta/2,delta/2,nbPoints); % disk
B{2} = shape.Ellipse(delta*1,delta/2,nbPoints); % ellipse
B{3} = shape.Flower(delta/2,delta/2,nbPoints); % flower
B{4} = shape.Imgshape([imagepath,'/A.png'], nbPoints); % A

% B{1} = shape.Ellipse(delta/2,delta/2,nbPoints); % disk
% B{2} = shape.Ellipse(delta*1,delta/2,nbPoints); % ellipse
% B{3} = shape.Flower(delta/2,delta/2,nbPoints); % flower
% B{4} = shape.Triangle(delta, pi/3, nbPoints); % triangle
% B{5} = shape.Rectangle(delta, delta, nbPoints); % square
% B{6} = shape.Rectangle(delta/2, delta, nbPoints); % rectangle
% B{7} = shape.Imgshape([imagepath,'/A.png'], nbPoints); % A
% B{8} = shape.Imgshape([imagepath,'/E.png'], nbPoints); % E


%% Compute multiscale time-dependent CGPT

disp('Computation of theoretical time dependent CGPTs...');

%%
% Parameters
ord = 2; % order of CGPT dictionary
scl = 2; % number of scales
Scl = 1.5.^(0:scl-1); % Scaling parameters
Ntime = 2^7; % time interval length, this seems not have influence on the result of matching
Tmax0 = 5; % Time duration of the pulse waveform at the scale 1
Tmax = zeros(1, scl); % Time duration at each scale
dt = zeros(1, scl);
Fmax = zeros(1, scl);
df = zeros(1, scl);
waveform = zeros(scl, Ntime);
freqform = zeros(scl, Ntime);

CGPTt = cell(length(B), scl);

for s = 1:scl
    % pulse waveform at the scale s    
    [waveform(s,:), dt(s), Tmax(s), freqform(s, :), df(s), Fmax(s)] = tools.make_pulse(Tmax0, Ntime, Scl(s));
end

for m=1:length(B)
    fprintf('Proceeding the shape %s...\n', B{m}.name_str);

    tic
    for s = 1:scl
        % pulse waveform at the scale s
        % [CGPTt{m,s}, dt, Tmax(s), ~] = asymp.CGPT.theoretical_CGPT_time_recon(B{m}, cnd, pmtt, ord, Tmax0, Ntime, Scl(s));
        
        % Equivalently, one can do the following (much slower but more accurate):
        [CGPTt{n,s}, dt, ~] = asymp.CGPT.theoretical_CGPT_time(B{n}, cnd, pmtt, ord, freqform(s,:), df(s), Tmax(s), Ntime); % dt depends only on df
    end
    toc
end

%%
% Compute multiscale shape descriptors

SDt = zeros(length(B), Ntime*scl); % Shape descriptor

for n=1:length(B)
    SDt(n,:) = dico.CGPT.ShapeDescriptor_PT_time(CGPTt(n,:), Scl, dt, Tmax(1));
    % SDt(n,:) = dico.CGPT.ShapeDescriptor_PT_time(CGPTt(n,:), Scl);
end

%% Show graphically the distinction of shape descriptors

%%
% Fix the disk as the reference shape
W = SDt(1,:); 

%%
% Compare each shape with the disk
figure;
plot(SDt(2,:)-W, 'r'); hold on;
plot(SDt(3,:)-W, 'g'); 
plot(SDt(4,:)-W, 'b'); 

%% Shape recognition

%%
% We generate first the unknown shape by transforming an element from the dictionary

n=2;
D = (B{n}<0.2)*2 + [1,1]';

fprintf('CGPTs of a transformed shape...\n');
CGPTtD = {};
tic
for s = 1:scl
    [CGPTtD{s}, ~] = asymp.CGPT.theoretical_CGPT_time_recon(D, cnd, pmtt, 2, Tmax0, Ntime, Scl(s));
end
toc

SDr = dico.CGPT.ShapeDescriptor_PT_time(CGPTtD, Scl);
norm(SDr-SDt(n,:))/norm(SDt(n,:))

figure; plot(SDr, 'r'); hold on; plot(SDt(n,:));
figure; plot(SDr-W, 'r'); hold on; plot(SDt(n,:)-W); 
% figure; plot(SDr-SDt(n,:));

[err, idx] = dico.CGPT.SD_Matching_time(SDr, SDt)
