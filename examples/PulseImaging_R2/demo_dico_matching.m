%% Matching in a dictionary
close all;
clear all;
addpath('~/SIES/');
matlabpool open;

%% Load the dictionary
aperture = 1/32;

pathname = '~/Data/';
% pathname = '/Volumes/ExFAT200G/Data/';
load([pathname, 'measurements/Pulse/Transformed/',num2str(aperture),'pi/data11_6scl.mat']);
% load([pathname, 'measurements/Pulse/Original/',num2str(aperture),'pi/data11_6scl.mat']);

%% Load the dictionary and construct shape descriptors
load([pathname,'/dico/Pulse/smalldico11_6scl.mat']);
Dico.names{7} = 'Rectangle 2'; Dico.names{11} = 'Ellipse 2';

% Bidx = 1:length(Dico.B); % index of shapes to be identified
% Bidx = [1, 2, 4:6, 8, 10:11]; % without flower, rectangle 2, L

% Bidx = [1:3, 5, 6, 9:11]; % without triangle, rectangle 2, A
  
% Bidx = [1:3, 5, 6, 9:11]; % without triangle, rectangle 2, A
Bidx = [1:3, 5, 6, 8, 9, 11]; % without triangle, rectangle 2, E

nbShapes = length(Bidx);
B = Dico.B(Bidx);
names = Dico.names(Bidx);

cnd = Dico.cnd(Bidx);
pmtt = Dico.pmtt(Bidx);

% index of scales
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

%% Dico-matching

Hrecon = @(data, nlvl)P.addnoise_recon(data, nlvl, 1, 10^5, 10^-8, 1, 'lsqr', op);

Hsd = @(CGPTt)dico.CGPT.ShapeDescriptor_PT_time(CGPTt, SD_method, Scl);

Hmatching = @(SDt, SDt_Dico)dico.CGPT.SD_Matching_time(SDt, SDt_Dico);

data = cell(nbShapes,1);
for m=1:nbShapes    
    data{m} = Data.data(Bidx(m), Sidx);
end

NLvls0 = 0.5:0.5:8;
NLvls = NLvls0(1:end);
nbNlv = length(NLvls);

nbExp = 250;
Mrate = cell(nbNlv, nbShapes);
Err0 = cell(nbNlv, nbShapes); 

parfor k = 1:nbNlv
    fprintf('Proceeding the noise level %f...\n', NLvls(k));
    [Mrate(k,:), Err0(k,:)] = dico.montecarlo_matching(data, SDt_Dico, nbExp, NLvls(k), Hrecon, ...
                                                      Hsd, Hmatching, 0);
end

%% Interpretation of the result We show in a bar figure the similarity between dictionary
% shape descriptors and the one reconstructed from data.

Err = zeros(nbNlv, nbShapes, nbShapes, nbScl); Idx = Err;

for k = 1:nbNlv
    for m = 1:nbShapes
        for s=1:nbScl
            Err(k, m, :, s) = Err0{k,m}(:,s);
            [~, Idx(k,m,:,s)] = sort(Err(k, m, :, s));
        end
    end
end

ss=4; kk=1; % choose the scale and noise level
Err = squeeze(Err(kk,:,:,ss)); 
Idx = squeeze(Idx(kk,:,:,ss)); 

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

%% Save the results

Res.Err = Err;
Res.Idx = Idx;
Res.Mrate = Mrate;
Res.NLvls = NLvls;
Res.nbExp = nbExp;

pathname = ['~/Data/outputs/Pulse/',num2str(aperture),'pi/'];
% pathname = ['/Volumes/ExFAT200G/Data/outputs/Pulse/',num2str(aperture),'pi/'];
mkdir(pathname);
fname = [pathname,'matching_data',num2str(length(B)),'_', num2str(nbScl),'scl_', num2str(NLvls(1)), 'nlvl.mat'];

save(fname,'Res','-v7.3');
fprintf('Data saved in %s\n', fname);

%% Plot Mrate
Rrate0 = zeros(nbNlv, nbShapes);
scl = 1;
NLvls1 = 0.5:0.25:8;
Rrate = zeros(length(NLvls1), nbShapes);

for n=1:nbShapes
	for k=1:nbNlv
		Rrate0(k,n) = Mrate{k,n}(scl);
	end
	
	Rrate(:,n) = interp1(NLvls, Rrate0(:,n), NLvls1, 'cubic'); % Interpolation if necessary
end
	
fig = figure; 
plot(NLvls1, Rrate(:,1), '-bo'); hold on;
plot(NLvls1, Rrate(:,2), '--ro');
plot(NLvls1, Rrate(:,3), '-g*');
plot(NLvls1, Rrate(:,4), '-bs');
plot(NLvls1, Rrate(:,5), '--rs');
plot(NLvls1, Rrate(:,6), '-g^');
plot(NLvls1, Rrate(:,7), '--bv');
plot(NLvls1, Rrate(:,8), '-.go');

legend(names);

ylim([-0.01, 1.01]);

%% old version
% Err0 = zeros(nbNlv, nbShapes, nbShapes, nbScl); Idx0 = Err0;

% % out = cell(nbShapes, nbExp, nbScl);

% % SDt = {}; SD = {};
% % out = cell(nbShapes, nbExp, nbScl);
% % data_noisy = cell(nbExp,1);
% % SDt = cell(nbExp,1);

% for k = 1:nbNlv
%     nlvl = NLvls(k); % noise level
%     fprintf('Proceeding the noise level %f...\n', nlvl);
    
%     for m = 1:nbShapes
%         tic

%         fprintf('Proceeding the shape %s...\n', names{m});

%         errtmp = zeros(nbExp, nbShapes, nbScl);
        
%         for p = 1:nbExp        
%             % fprintf('...Proceeding the %d-th trial of %d...\n', p, nbExp);
            
%             CGPTt = {};

%             for s = 1:nbScl
%                 ss = Sidx(s);
                
%                 data_noisy = P.add_white_noise_global(data{m,ss}, nlvl);
%                 out = P.reconstruct_CGPT(data_noisy.MSR_noisy, 1, 10^5, 10^-8, 1, 'lsqr', op);
%                 % out = P.reconstruct_CGPT(data_noisy.MSR_noisy, 1); % pinv
                
%                 CGPTt{s} = tools.cell2mat3D(out.CGPT);
%             end

%             % [SDt{m,p}, SD{m,p}] = dico.CGPT.ShapeDescriptor_PT_time(CGPTt(m,p,:), SD_method, Scl);
%             % SDt = SDt(Dico.extrema, :);
%             % [errtmp(p,:,:), ~] = dico.CGPT.SD_Matching_time(SDt{m,p}, SDt_Dico);

%             [SDt, ~] = dico.CGPT.ShapeDescriptor_PT_time(CGPTt, SD_method, Scl);
%             [errtmp(p,:,:), ~] = dico.CGPT.SD_Matching_time(SDt, SDt_Dico);
%         end
        
%         Err0(k,m,:,:) = mean(errtmp,1);
%         for s=1:nbScl
%             [~, Idx0(k,m,:,s)] = sort(Err0(k,m,:,s));
%         end
%         toc
%     end
% end

% Res.Err0 = Err0;
% Res.Idx0 = Idx0;
% Res.nbExp = nbExp;
% Res.NLvls = NLvls;

% % pathname = ['~/Data/measurements/Pulse/Transformed/',num2str(aperture),'pi/'];
% pathname = ['~/Data/outputs/Pulse/',num2str(aperture),'pi/'];
% % pathname = ['/Volumes/ExFAT200G/Data/measurements/Pulse/Transformed/',num2str(aperture),'pi/'];
% % pathname = ['/Volumes/ExFAT200G/Data/measurements/Pulse/Original/',num2str(aperture),'pi/'];
% mkdir(pathname);
% fname = [pathname,'data',num2str(length(B)),'_', num2str(nbScl),'scl_', num2str(NLvls(1)), 'nlvl.mat'];

% save(fname,'Res','-v7.3');
% fprintf('Data saved in %s\n', fname);

% cc
