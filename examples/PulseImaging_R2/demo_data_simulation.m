%% Data simulation for Pulse imaging.
% For shapes of the dictionary, we apply some perturbation and simulate
% the MSR stream data. The parameter sets here need to be the same as in
% make_dico.m.

%%
clear all;
close all;
addpath('~/SIES');
% matlabpool;

%% Load the dictionary and construct shape descriptors
% pathname = '/Volumes/ExFAT200G/Data/';
pathname = '~/Data/';
dico_name = [pathname,'/dico/Pulse/smalldico11_6scl.mat'];
load(dico_name);

Bidx = 1:length(Dico.B); % index of shapes that data will be simulated
B = Dico.B(Bidx);
cnd = Dico.cnd(Bidx);
pmtt = Dico.pmtt(Bidx);

%% Parameters and waveforms

Sidx = 1:length(Dico.Scl); % index of scales that data will be simulated
Scl = Dico.Scl(Sidx);
nbScl = length(Scl); % number of scales

Ntime = 2^9; % time interval length
Tmax = zeros(1, nbScl);
dt = zeros(1, nbScl);
waveform = zeros(nbScl, Ntime);

for s = 1:nbScl
    % pulse waveform at the scale s
    [waveform(s,:), dt(s), Tmax(s), ~] = tools.make_pulse(Dico.Tmax0, Ntime, Scl(s));
end

%% Data simulation

Ns = 50; % Number of sources

rot = pi/3; sca = 1.5; trl = 0.1*[1, 1]';
% rot = 0; sca = 1; trl = 0.*[1, 1]';

mradius = max(12*(sca/2+norm(trl)), 10);

Aperture = [1/32, 1/16, 1/8, 1/4, 1/2, 1, 2];

for aa=1:length(Aperture)
    aperture = Aperture(aa);
    
    fprintf('Aperture angle: %f\n', aperture);
    
    cfg = acq.Coincided([0,0]', mradius, Ns, [1, aperture*pi, 2*pi], false, [1,-1], 0.01);
    
    data = cell(length(B), nbScl);
    
    % Compute time-dependent CGPT
    for s = 1:nbScl
        fprintf('...Proceeding the scale %f...\n', Scl(s));
        
        parfor n=1:length(B) % iteration on the shape
            fprintf('......Proceeding the shape %s...\n', B{n}.name_str);
        
            D = (B{n}<rot)*sca + trl;
            
            P = PDE.PulseImaging_R2(D, cnd(n), pmtt(n), waveform(s,:), dt(s), cfg);
            % figure; plot(P); axis image;
            
            data{n,s} = P.data_simulation();
        end
    end
    
    Data = [];
    Data.data = data;
    Data.B = B;
    Data.cnd = cnd;
    Data.pmtt = pmtt;
    
    Data.cfg = cfg; % We keep only the configuration
    Data.Scl = Scl;
    Data.Sidx = Sidx;
    Data.dico_name = dico_name;

    Data.Ntime = Ntime;
    Data.Tmax = Tmax;
    Data.dt = dt;
    Data.waveform = waveform;
    
    Data.rot = rot;
    Data.sca = sca;
    Data.trl = trl;
    
    pathname = ['~/Data/measurements/Pulse/Transformed/',num2str(aperture),'pi/'];
    % pathname = ['/Volumes/ExFAT200G/Data/measurements/Pulse/Transformed/',num2str(aperture),'pi/'];
    % pathname = ['/Volumes/ExFAT200G/Data/measurements/Pulse/Original/',num2str(aperture),'pi/'];
    mkdir(pathname);
    fname = [pathname,'data',num2str(length(B)),'_', num2str(nbScl),'scl.mat'];
    
    save(fname,'Data','-v7.3');
    fprintf('Data saved in %s\n', fname);
end
