%% Matching in a dictionary
close all;
clear all;
addpath('~/SIES/');

%% Load the dictionary
aperture = 1;
% pathname = '/Volumes/Yue/Data/measurements/Pulse/Transformed/0.03125pi/';
% pathname = '~/Data/measurements/Pulse/Transformed/1pi/';

pathname = '/Volumes/ExFAT200G/Data/';
% load([pathname, 'measurements/Pulse/Transformed/',num2str(aperture),'pi/data9_4scl.mat']);
load([pathname, 'measurements/Pulse/Original/',num2str(aperture),'pi/data11_6scl.mat']);

%% Load the dictionary and construct shape descriptors
load([pathname,'/dico/Pulse/smalldico11_6scl.mat']);

% Bidx = 1:length(Dico.B); % index of shapes to be identified
Bidx = [1:2, 4,5, 7:10];
% Bidx = [1, 2, 3, 9];
nbShapes = length(Bidx);
B = Dico.B(Bidx);
names = Dico.names(Bidx);

cnd = Dico.cnd(Bidx);
pmtt = Dico.pmtt(Bidx);

Sidx =  1:length(Dico.Scl); % index of scales
% Sidx =  11:15; % index of scales
Scl = Dico.Scl(Sidx);
nbScl = length(Scl);

SDt_Dico = {}; % Shape descriptor
SD_Dico = {}; % Shape descriptor

SD_method = 2; % Method for construction of shape descriptors

for m=1:nbShapes % iteration on the shape
    for s = 1:nbScl
        ss = Sidx(s);
        [CGPTt0{m,s}, dt(s)] = asymp.CGPT.CGPT_time_truncation(Dico.CGPTt0{Bidx(m),ss}, Dico.dt0(ss), Data.Tmax(s), Data.Ntime);
    end

    [SDt_Dico{m}, SD_Dico{m}] = dico.CGPT.ShapeDescriptor_PT_time(CGPTt0(m,:), SD_method, Scl);

    % % Or use only the extrema as shape descriptors (less robust to noise)
    %     toto = dico.CGPT.ShapeDescriptor_PT_time(CGPTt, Scl, SD_method);
    %     SDt_Dico{m} = toto(Dico.extrema{1}, :);

end

%% Test the dictionary of shape descriptors
Err = zeros(nbShapes, nbShapes); Idx = Err;

for m = 1:nbShapes
    [Err(m,:), Idx(m,:)] = dico.CGPT.SD_Matching_time(SDt_Dico{m}, SDt_Dico);
    % [Err(m,:), Idx(m,:)] = dico.CGPT.SD_Matching_time(SDt(m,1), SDt(:,1));
end

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
fig=figure; plot(P); axis image; xlim([-1,1]*12); ylim([-1,1]*12);
% fig=figure; plot(P); axis image; xlim([-1,1]*8); ylim([-0.2,1]*8);
% saveas(fig, '../figures/limview1pi.eps', 'psc2');

%% Dico-matching
nlvl = 3; % noise level
nbExp = 1;

data = Data.data(Bidx, :);

Err = zeros(nbShapes, nbShapes); Idx = Err;

out = cell(nbShapes, nbExp, nbScl);
CGPTt = cell(nbShapes, nbExp, nbScl);
SDt = {}; SD = {};
% out = cell(nbShapes, 1);
% data_noisy = cell(nbExp,1);
% SDt = cell(nbExp,1);

for m = 1:nbShapes
    fprintf('Proceeding the shape %s...\n', names{m});

    errtmp = zeros(nbExp, nbShapes);

    for i = 1:nbExp
        for s = 1:nbScl
            fprintf('...Proceeding the scale %f...\n', Scl(s));
            data_noisy = P.add_white_noise_global(data{m,s}, nlvl);
            % out{m,i,s} = P.reconstruct_CGPT(data_noisy.MSR_noisy, 1, 10^5, 10^-8, 1, 'lsqr');
            out{m,i,s} = P.reconstruct_CGPT(data_noisy.MSR_noisy, 1); % pinv
            
            CGPTt{m,i,s} = tools.cell2mat3D(out{m,i,s}.CGPT);
        end
        
        [SDt{m,i}, SD{m,i}] = dico.CGPT.ShapeDescriptor_PT_time(CGPTt(m,i,:), SD_method, Scl);
        % SDt = SDt(Dico.extrema, :);
        [errtmp(i,:), ~] = dico.CGPT.SD_Matching_time(SDt{m,i}, SDt_Dico, 1:nbScl);
    end
    
    Err(m,:) = mean(errtmp,1);
    [~, Idx(m,:)] = sort(Err(m,:));
end

% Interpretation of the result We show in a bar figure the similarity between dictionary
% shape descriptors and the one reconstructed from data.

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

%% Compare the MSR
n=8;
s=4; 

data_noisy = P.add_white_noise_global(data{n,s}, nlvl);
MSR0 = tools.cell2mat3D(data_noisy.MSR);
MSR1 = tools.cell2mat3D(data_noisy.MSR_noisy);

rr=1; cc= 1;
toto0 = MSR0(rr,cc,:); toto1 = MSR1(rr,cc,:);
xt = (0:Data.Ntime-1)*Data.dt(s);
fig = figure; plot(xt, squeeze(toto0), 'LineWidth',2); hold on; plot(xt, squeeze(toto1), 'r');
xlabel('Time'); %ylim(1.5*[-1,1]*1e-3);

% saveas(fig, '~/Writings/PulseImaging/figures/MSR_noisy_scl1_nlvl0p5.eps' ,'psc2');

%% Error of the reconstruction
toto1 = zeros(1, Data.Ntime);
toto2 = zeros(1, Data.Ntime);

n=3; s=1;

for i = 1:nbExp
    for t=1:Data.Ntime
        toto1(t) = toto1(t) + out{n,i,s}.res{t};
        toto2(t) = toto2(t) + out{n,i,s}.rres{t};
    end
end

aerr = toto1 / nbExp;
rerr = toto2 / nbExp;
figure; plot(aerr);
figure; plot(rerr);

%% Compare the CGPT (meaningful only if the object is identical to some shape of the dictionary, without scaling or rotation)
% 
n=8; s=4;
row = 2; col = 2;
figure; plot(squeeze(CGPTt0{n,s}(row, col,:))); hold on; plot(squeeze(CGPTt{n,1,s}(row, col, :)), 'r')
toto = CGPTt{n,1,s}(row, col,:) - CGPTt0{n,s}(row, col,:);
figure; plot(squeeze(toto(1,1,:)));

%%
n=8;
% figure; plot(SDt{8,1}(:,s)); hold on; plot(SDt_Dico{8}(:,s), 'r')
t1 = squeeze(SD{n,1}(:,s,1)); s1 = squeeze(SD_Dico{n}(:,s,1));
t2 = squeeze(SD{n,1}(:,s,2)); s2 = squeeze(SD_Dico{n}(:,s,2));
figure; plot(t1); hold on; plot(s1, 'r')
figure; plot(t2); hold on; plot(s2, 'r')
% figure; plot(t1.^2+t2.^2); hold on; plot(s1.^2+s2.^2, 'r')
% toto = SDt{n,1}(:,s) - SDt_Dico{n}(:,s);
% figure; plot(toto);

%%
for n=1:nbShapes
    err = 0;
    for i=1:nbExp
        toto = CGPTt{n,i,s}(1,1,:) - CGPTt0{n,s}(1,1,:);
        err = err + norm(squeeze(toto(1,1,:)), 'fro');
    end
    err = err/nbExp;
    err/norm(squeeze(CGPTt0{n,s}(1,1,:)),'fro')
end
