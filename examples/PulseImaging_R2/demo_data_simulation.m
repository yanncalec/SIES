%% Data simulation for Pulse imaging.
% For each shape of the dictionary, we apply some perturbation and simulate
% the MSR stream data. The parameter sets here need to be compatible with
% the ones in make_dico.m

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
cnd = 2;
pmtt = 2;

%% Make dictionary

disp('Construction of the dictionary...');

% B{1} = shape.Flower(delta/2,delta/2,nbPoints); % flower

% B{1} = shape.Ellipse(delta/2,delta/2,nbPoints); % disk
% B{2} = shape.Ellipse(delta*1,delta/2,nbPoints); % ellipse
% B{3} = shape.Flower(delta/2,delta/2,nbPoints); % flower
% B{4} = shape.Imgshape([imagepath,'/A.png'], nbPoints); % A

B{1} = shape.Ellipse(delta/2,delta/2,nbPoints); % disk
B{2} = shape.Ellipse(delta*1,delta/2,nbPoints); % ellipse
B{3} = shape.Flower(delta/2,delta/2,nbPoints); % flower
B{4} = shape.Triangle(delta, pi/3, nbPoints); % triangle
B{5} = shape.Rectangle(delta, delta, nbPoints); % square
B{6} = shape.Rectangle(delta/2, delta, nbPoints); % rectangle
B{7} = shape.Imgshape([imagepath,'/A.png'], nbPoints); % A
B{8} = shape.Imgshape([imagepath,'/E.png'], nbPoints); % E

%% Compute multiscale time-dependent CGPT

disp('Computation of theoretical time dependent CGPTs...');
%%
% Parameters
Scl = 1.5.^[0];
scl = length(Scl); % number of scales

Ntime = 2^9; % time interval length
Tmax0 = 5;
Tmax = zeros(1, scl);
dt = zeros(1, scl);
waveform = zeros(scl, Ntime);

for s = 1:scl
    % pulse waveform at the scale s    
    [waveform(s,:), dt(s), Tmax(s), ~] = tools.make_pulse(Tmax0, Ntime, Scl(s));
end

%%
% Data simulation

% rot = pi/3; sca = 1.5; trl = [0.5, 0.5]';
% tflag = true;
% mradius = 5*(sca/2+norm(trl)); 

rot = 0; sca = 1.; trl = [0., 0.]';
tflag = false;
mradius = 3.5;

% Aperture = [0.125, 0.25, 0.5, 0.75, 1, 2];
Aperture = [1];

for aa=1:length(Aperture)
    aperture = Aperture(aa);
	fprintf('Aperture angle: %f\n', aperture);

	cfg = acq.Coincided([0,0]', mradius, 50, [1, aperture*pi, 2*pi], false, [1,-1]);
    
    data = cell(length(B), scl);

    % Compute time-dependent CGPT    
    tic
    for n=1:length(B) % iteration on the shape
        fprintf('Proceeding the shape %s...\n', B{n}.name_str);               
        
        D = (B{n}<rot)*sca + trl;

        for s = 1:scl
            P = PDE.PulseImaging_R2(D, cnd, pmtt, waveform(s,:), dt(s), cfg);
            % figure; plot(P); axis image;
            
            data{n, s} = P.data_simulation();
        end
    end
    
    Data = [];
    Data.data = data;
    Data.B = B;
    Data.cfg = cfg;

    Data.Scl = Scl;
    Data.Tmax0 = Tmax0;
    Data.Ntime = Ntime;
    Data.Tmax = Tmax;
    Data.dt = dt;
    Data.waveform = waveform;
    
    Data.rot = rot;
    Data.sca = sca;
    Data.trl = trl;
    
    if ~tflag
        fname = ['/Volumes/Yue/Pulse/Original/'
            /Data/measurements/Pulse/Original/',num2str(aperture),'pi/data',num2str(length(B)),'_', num2str(scl),'scl.mat'];
    else
        fname = ['~/Data/measurements/Pulse/Transformed/',num2str(aperture),'pi/data',num2str(length(B)),'_', num2str(scl),'scl.mat'];
    end
    
    save(fname,'Data','-v7.3');
    fprintf('Data saved in %s\n', fname);
    toc
end

% fig = figure; plot(P); axis image;
% xlim(1.1*[-mradius,mradius]); ylim(1.1*[-mradius,mradius]);
