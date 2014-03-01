%% Data simulation for all shapes in a small dictionary using multiple frequencies

%% 
% Load the dictionary
load ~/Data/dico/SCT/smalldico.mat;

Nf = length(Dico.freq);
sclrange = Dico.sclrange;
sfrange = Dico.sfrange;

pmtt = Dico.pmtt;
pmeb = Dico.pmeb;
pmtt_bg = Dico.pmtt_bg;
pmeb_bg = Dico.pmeb_bg;

%%
% Number of sources (and receivers)
N0 = 91; 

%%
% Transform parameters of a standard shape
scl = 1.5; trl = [-0.5,0.5]'; rtn = pi/3;
% scl = 1.75; trl = [-0.5,0.5]'; rtn = pi/3;

%%
% Scanning frequencies
Nsf = floor(Nf/sclrange(2)); % number of scanning freqs
sfreq = linspace(sfrange(1), sfrange(2), Nsf);

%%
% Make the acquisition configuration. With the class |acq.Planewave|, we make equally distributed
% sources/receivers on a circle of center at [0,0]' with radius |mradius|.

N0 = 91; % Number of sources (and receivers)
cfg = acq.Planewave([0,0]', 3, N0, N0, [1, 2*pi, 2*pi]);

disp('Data simulation...');

for n = 1:length(Dico.B)
    fprintf('Processing the shape: %s...\n', Dico.B{n}.name_str);

    D = (Dico.B{n}<rtn) * scl + trl;
    
    P = PDE.Helmholtz_R2(D, pmtt, pmeb, pmtt_bg, pmeb_bg, cfg); 
    % figure; plot(P, 'LineWidth', 1); axis image;

    out{n} = P.data_simulation(sfreq);    
end

Data.out = out;
Data.scl = scl;
Data.trl = trl;
Data.rtn = rtn;
Data.sfreq = sfreq;
Data.cfg = cfg;

save('~/Data/out/SCT/MSR.mat', 'Data', '-v7.3');
