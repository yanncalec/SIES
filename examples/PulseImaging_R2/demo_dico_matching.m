%% Monte-carlo test of dictionary matching algorithm
% This script takes as in put the dictionary of time dependent CGPTs computed by
% demo_make_dico.m, and the MSR stream simulated by demo_data_simulation.m,
% reconstruct the CGPTs from data and apply the matching algorithm at
% various level of noise.
% 

close all;
clear all;
addpath('~/SIES/');

%% Load the simulated MSR data
aperture = 1/16;

pathname = '/Volumes/Yue_Fat32/Data/';
load([pathname, 'measurements/Pulse/Transformed/',num2str(aperture),'pi/data11_6scl.mat']);

%% Load the dictionary and construct shape descriptors
load([pathname,'/dico/Pulse/smalldico11_6scl.mat']);
Dico.names{7} = 'Rectangle 2'; Dico.names{11} = 'Ellipse 2';

% Choose a subset of shapes that will be used for matching
Bidx = [1:3, 5, 6, 8, 9, 11]; % without triangle, rectangle 2, E
% Bidx = 1:length(Dico.B); % index of shapes to be identified
% Bidx = [1, 2, 4:6, 8, 10:11]; % without flower, rectangle 2, L
% Bidx = [1:3, 5, 6, 9:11]; % without triangle, rectangle 2, A
% Bidx = [1:3, 5, 6, 9:11]; % without triangle, rectangle 2, A

nbShapes = length(Bidx);
B = Dico.B(Bidx);
names = Dico.names(Bidx);

cnd = Dico.cnd(Bidx);
pmtt = Dico.pmtt(Bidx);

% Choose a subset of scales that will be used for matching
% Sidx =  1:length(Dico.Scl); 
Sidx =  2:5;
Scl = Dico.Scl(Sidx);
nbScl = length(Scl);

SDt_Dico = {}; % Shape descriptor

% SD_Dico = {}; 

SD_method = 2; % Method for construction of shape descriptors

CGPTt0 = {};

for m=1:nbShapes % iteration on the shape
    for s = 1:nbScl
        ss = Sidx(s);
        [CGPTt0{m,s}, dt(s)] = asymp.CGPT.CGPT_time_truncation(Dico.CGPTt0{Bidx(m),ss}, Dico.dt0(ss), Data.Tmax(ss), Data.Ntime);
    end

    % [SDt_Dico{m}, SD_Dico{m}] = dico.CGPT.ShapeDescriptor_PT_time(CGPTt0(m,:), SD_method, Scl);
    [SDt_Dico{m}, ~] = dico.CGPT.ShapeDescriptor_PT_time(CGPTt0(m,:), SD_method, Scl);

    % % Or use only the extrema as shape descriptors (less robust to noise)
    %     toto = dico.CGPT.ShapeDescriptor_PT_time(CGPTt, Scl, SD_method);
    %     SDt_Dico{m} = toto(Dico.extrema{1}, :);

end

%% Test the dictionary of shape descriptors
Err0 = zeros(nbShapes, nbShapes, nbScl); Idx0 = Err0;

for m = 1:nbShapes
    [Err0(m,:,:), ~, Idx0(m,:,:)] = dico.CGPT.SD_Matching_time(SDt_Dico{m}, SDt_Dico);
end

Err = Err0(:,:,end); Idx = Idx0(:,:,end); % use all scales

fig=figure;
bar(Err, 'facecolor', 'none'); 
set(gca, 'XTickLabel', names, 'XTick',1:nbShapes);
hold on; 

% Identification result is plotted in green, if the result is wrong, the
% true shape is marked in red.
bar(eye(nbShapes).*Err, 'r'); 

toto=eye(nbShapes); 
idx = Idx(:,1); 
bar(toto(idx, :).*Err, 'g'); 

%% Parameters 

% Construct the PDE environment. Only the geometrical settings cfg in P
% will be used in the reconstruction.
cfg = Data.cfg;
P = PDE.PulseImaging_R2((B{1}<pi/3)*1.5+0.1*[1,1]', cnd(2), pmtt(2), Data.waveform(2,:), Data.dt(2), cfg); 
op = PDE.Conductivity_R2.make_linop_CGPT(cfg, 1, 1); % Construct the linear operator L

% fig=figure; plot(P); axis image; xlim([-1,1]*12); ylim([-1,1]*12);
% fig=figure; plot(P); axis image; xlim([-1,1]*8); ylim([-0.2,1]*8);
% saveas(fig, '../figures/limview1pi.eps', 'psc2');

%% Monte-carlo test of dictionary matching

Hrecon = @(data, nlvl)P.addnoise_recon(data, nlvl, 1, 10^5, 10^-8, 1, 'lsqr', op);

Hsd = @(CGPTt)dico.CGPT.ShapeDescriptor_PT_time(CGPTt, SD_method, Scl);

Hmatching = @(SDt, SDt_Dico)dico.CGPT.SD_Matching_time(SDt, SDt_Dico);

data = cell(nbShapes,1);
for m=1:nbShapes    
    data{m} = Data.data(Bidx(m), Sidx);
end

NLvls0 = [2]; 
NLvls = NLvls0(1:end); % level (percentage) of noise
nbNlv = length(NLvls);

nbExp = 1;
Mrate = cell(nbNlv, nbShapes);
Err0 = cell(nbNlv, nbShapes); 

for k = 1:nbNlv
    fprintf('Proceeding the noise level %f...\n', NLvls(k));
    [Mrate(k,:), Err0(k,:)] = dico.montecarlo_matching(data, SDt_Dico, nbExp, NLvls(k), Hrecon, ...
                                                      Hsd, Hmatching, 0);
end

%% Interpretation of the result We show in a bar figure the similarity between dictionary
% shape descriptors and the one reconstructed from data.

Err_all = zeros(nbNlv, nbShapes, nbShapes, nbScl); Idx_all = Err_all;

for k = 1:nbNlv
    for m = 1:nbShapes
        for s=1:nbScl
            Err_all(k, m, :, s) = Err0{k,m}(:,s);
            [~, Idx_all(k,m,:,s)] = sort(Err_all(k, m, :, s));
        end
    end
end

ss=4; kk=1; % choose the scale and noise level
Err = squeeze(Err_all(kk,:,:,ss)); 
Idx = squeeze(Idx_all(kk,:,:,ss)); 

fig1= figure; 
bar(Err, 'facecolor', 'none'); 
set(gca, 'XTickLabel', names, 'XTick',1:nbShapes);
hold on; 

% Identification result is plotted in green, if the result is wrong, the
% true shape is marked in red.
bar(eye(nbShapes).*Err, 'r'); 

toto = eye(nbShapes); 
idx = Idx(:,1); 
bar(toto(idx, :).*Err, 'g'); 
ylim([0, 0.375]);
saveas(fig1, '../figures/ErrorBar_pi16_nlvl2.eps', 'psc2');

%% Save the results

Res.Err_all = Err_all;
Res.Idx_all = Idx_all;
Res.Mrate = Mrate;
Res.NLvls = NLvls;
Res.nbExp = nbExp;

pathname = ['~/Data/outputs/Pulse/',num2str(aperture),'pi/'];
% pathname = ['/Volumes/ExFAT200G/Data/outputs/Pulse/',num2str(aperture),'pi/'];
mkdir(pathname);
fname = [pathname,'matching_data',num2str(length(B)),'_', num2str(nbScl),'scl_', num2str(NLvls(1)), 'nlvl.mat'];

save(fname,'Res','-v7.3');
fprintf('Data saved in %s\n', fname);

%% Plot results of monte-carlo test
MC0 = zeros(nbNlv, nbShapes);
scl = 1;
NLvls1 = 0.5:0.25:8;
MC = zeros(length(NLvls1), nbShapes);

for n=1:nbShapes
	for k=1:nbNlv
		MC0(k,n) = Mrate{k,n}(scl);
	end
	
	MC(:,n) = interp1(NLvls, MC0(:,n), NLvls1, 'cubic'); % Interpolation if necessary
end
	
fig = figure; 
plot(NLvls1, MC(:,1), '-bo'); hold on;
plot(NLvls1, MC(:,2), '--ro');
plot(NLvls1, MC(:,3), '-g*');
plot(NLvls1, MC(:,4), '-bs');
plot(NLvls1, MC(:,5), '--rs');
plot(NLvls1, MC(:,6), '-g^');
plot(NLvls1, MC(:,7), '--bv');
plot(NLvls1, MC(:,8), '-.go');

legend(names);

ylim([-0.01, 1.01]);
