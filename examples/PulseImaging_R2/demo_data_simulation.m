%% Data simulation for Pulse imaging.
% For shapes of the dictionary, we apply some perturbation and simulate
% the MSR stream data. The parameter sets here need to be the same as in
% make_dico.m.

%%
clear all;
close all;
addpath('~/SIES');

%% Make dictionary
load ~/Data/dico/Pulse/smalldico9_6scl.mat;

Bidx = 1:length(Dico.B); % index of shapes that data will be simulated
B = Dico.B(Bidx);
cnd = Dico.cnd(Bidx);
pmtt = Dico.pmtt(Bidx);

%% Compute multiscale time-dependent CGPT

disp('Computation of theoretical time dependent CGPTs...');
%%
% Parameters
scl = length(Dico.Scl); % number of scales

Ntime = 2^9; % time interval length
Tmax = zeros(1, scl);
dt = zeros(1, scl);
waveform = zeros(scl, Ntime);

for s = 1:scl
    % pulse waveform at the scale s    
    [waveform(s,:), dt(s), Tmax(s), ~] = tools.make_pulse(Dico.Tmax0, Ntime, Dico.Scl(s));
end

%%
% Data simulation

Ns = 50; % Number of sources

rot = pi/3; sca = 1.5; trl = [0.5, 0.5]';
mradius = 5*(sca/2+norm(trl)); 

Aperture = [0.125, 0.25, 0.5, 0.75, 1, 2];

for aa=1:length(Aperture)
    aperture = Aperture(aa);
	
    fprintf('Aperture angle: %f\n', aperture);

	cfg = acq.Coincided([0,0]', mradius, Ns, [1, aperture*pi, 2*pi], false, [1,-1]);
    
    data = cell(length(B), scl);

    % Compute time-dependent CGPT    
    tic
    for n=1:length(B) % iteration on the shape
        fprintf('Proceeding the shape %s...\n', B{n}.name_str);
        
        D = (B{n}<rot)*sca + trl;

        for s = 1:scl
            P = PDE.PulseImaging_R2(D, cnd(n), pmtt(n), waveform(s,:), dt(s), cfg);
            % figure; plot(P); axis image;
            
            data{n, s} = P.data_simulation();
        end
    end
    
    Data = [];
    Data.data = data;
    Data.B = B;
    Data.cnd = cnd;
    Data.pmtt = pmtt;
    
    Data.cfg = cfg; % We keep only the configuration
    Data.Scl = Dico.Scl;
    Data.Ntime = Ntime;
    Data.Tmax = Tmax;
    Data.dt = dt;
    Data.waveform = waveform;
    
    Data.rot = rot;
    Data.sca = sca;
    Data.trl = trl;
    
    fname = ['~/Data/measurements/Pulse/Transformed/',num2str(aperture),'pi/data',num2str(length(B)),'_', num2str(scl),'scl.mat'];
    
    save(fname,'Data','-v7.3');
    fprintf('Data saved in %s\n', fname);
    toc
end

%% Some manip
Aperture = [0.125, 0.25, 0.5, 0.75, 1, 2];

for aa=1:length(Aperture)
    aperture = Aperture(aa);
    pathname = ['/Volumes/Yue_Fat32/Pulse/Transformed/', num2str(aperture),'pi/'];
    % pathname = '/Volumes/Yue_Fat32/Pulse/Transformed/2pi/';
    fname = [pathname, 'data1_6scl.mat'];    
    load(fname);    
    Data1 = Data;
    
    fname = [pathname, 'data8_6scl.mat'];    
    load(fname);

    Data.data(9,:) = Data1.data;
    Data.B = Dico.B;
    Data.cnd = Dico.cnd;
    Data.pmtt = Dico.pmtt;
    
    Data.Scl = Dico.Scl;
    fname = [pathname,'data',num2str(length(Dico.B)),'_', num2str(scl),'scl.mat'];
    
    save(fname,'Data','-v7.3');
    fprintf('Data saved in %s\n', fname);    
end