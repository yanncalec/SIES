%% Matching in a dictionary
close all;
clear all;
addpath('~/SIES/');

%% Load the dictionary

% pathname = '/Volumes/Yue/Data/measurements/Pulse/Transformed/0.03125pi/';
pathname = '/Volumes/Yue/Data/measurements/Pulse/Transformed/0.0625pi/';
% pathname = '/Volumes/Yue/Data/measurements/Pulse/Transformed/1pi/';
load([pathname, 'data9_6scl.mat']);

%% Load the dictionary and construct shape descriptors
% load ~/Data/dico/Pulse/smalldico9_6scl.mat;

Bidx = 1:length(Dico.B); % index of shapes to be identified
% Bidx = [1:2, 4:9];
nbShapes = length(Bidx);
B = Dico.B(Bidx);
names = Dico.names(Bidx);

cnd = Dico.cnd(Bidx);
pmtt = Dico.pmtt(Bidx);

Sidx =  length(Dico.Scl):-1:1; % index of scales
nbScl = length(Sidx);
Scl = Dico.Scl(Sidx);

SDt_Dico = {}; % Shape descriptor

SD_method = 2; % Method for construction of shape descriptors

for m=1:nbShapes % iteration on the shape
    for s = 1:length(Scl)
        % [CGPTt0{m,s}, dt(s)] = asymp.CGPT.CGPT_time_truncation(Dico.CGPTt0{Bidx(m),Sidx(s)}, Dico.dt0(s), Data.Tmax(s), Data.Ntime);
        [CGPTt0{m,s}, dt(s)] = asymp.CGPT.CGPT_time_truncation(Dico.CGPTt0{Bidx(m),Sidx(s)}, Dico.dt0(s), Dico.Tmax(s), Dico.Ntime);
    end

    SDt_Dico{m} = dico.CGPT.ShapeDescriptor_PT_time(CGPTt0(m,:), Scl, SD_method);

    % % Or use only the extrema as shape descriptors (less robust to noise)
    %     toto = dico.CGPT.ShapeDescriptor_PT_time(CGPTt, Scl, SD_method);
    %     SDt_Dico{m} = toto(Dico.extrema{1}, :);

end

%% Parameters 

% Construct the PDE environment. Only the geometrical settings cfg in P
% will be used in the reconstruction.
cfg = Data.cfg;
P = PDE.PulseImaging_R2((B{6}<pi/3)*1.5+[0.5,0.5]', cnd(2), pmtt(2), Data.waveform(2,:), Data.dt(2), cfg); 
% fig=figure; plot(P); axis image; xlim([-1,1]*8); ylim([-1,1]*8);
fig=figure; plot(P); axis image; xlim([-1,1]*8); ylim([-0.2,1]*8);
% saveas(fig, '../figures/limview1pi.eps', 'psc2');

%% Dico-matching
nlvl = 1; % noise level
nbExp = 1;

data = Data.data(Bidx, Sidx);

Err = zeros(nbShapes, nbShapes); Idx = Err;

CGPTt = cell(nbShapes, nbScl);
out = cell(nbShapes, nbExp, nbScl);

for m = 1:nbShapes
    errtmp = zeros(nbExp, nbShapes);
    %idx = zeros(nbExp, nbShapes);
    
    for i = 1:nbExp
        for s = 1:nbScl
            data_noisy = P.add_white_noise_global(data{m,s}, nlvl);
            out{m,i,s} = P.reconstruct_CGPT(data_noisy.MSR_noisy, 1, 10^5, 10^-5, 1, 'lsqr');
            % out{m,i,s} = P.reconstruct_CGPT(data_noisy.MSR_noisy, 1); % pinv
            
            CGPTt{m,s} = tools.cell2mat3D(out{m,i,s}.CGPT);
        end
        
        SDt = dico.CGPT.ShapeDescriptor_PT_time(CGPTt(m,:), Scl, SD_method);
        % SDt = SDt(Dico.extrema, :);
        [errtmp(i,:), ~] = dico.CGPT.SD_Matching_time(SDt, SDt_Dico, 1:6);
    end
    
    Err(m,:) = mean(errtmp,1);
    [~, Idx(m,:)] = sort(Err(m,:));
end

% Interpretation of the result
% We show in a bar figure the similarity between dictionary shape descriptors and the one reconstructed from data.

fig1= figure; 
bar(Err, 'facecolor', 'none'); 
set(gca, 'XTickLabel', names, 'XTick',1:nbShapes);
hold on; 

% Identification result is plotted in green, if the result is wrong, the
% true shape is marked in red.
bar(eye(nbShapes).*Err, 'r'); 

toto=eye(nbShapes); 
idx = Idx(:,1); 
bar(toto(idx, :).*Err, 'g'); 
cc
%% Compare the MSR
n=2;
s=1; 

data_noisy = P.add_white_noise_global(data{n,s}, nlvl);
MSR0 = tools.cell2mat3D(data_noisy.MSR);
MSR1 = tools.cell2mat3D(data_noisy.MSR_noisy);

rr=1; cc= 1;
toto0 = MSR0(rr,cc,:); toto1 = MSR1(rr,cc,:);
xt = (0:Data.Ntime-1)*Data.dt(s);
fig = figure; plot(xt, squeeze(toto0), 'LineWidth',2); hold on; plot(xt, squeeze(toto1), 'r');
xlabel('Time'); ylim(1.5*[-1,1]*1e-3);

% saveas(fig, '~/Writings/PulseImaging/figures/MSR_noisy_scl1_nlvl0p5.eps' ,'psc2');

%% Error of the reconstruction
toto1 = zeros(1, Data.Ntime);
toto2 = zeros(1, Data.Ntime);

for i = 1:nbExp
    for t=1:Data.Ntime
        toto1(t) = toto1(t) + out{n,i,s}.res{t};
        toto2(t) = toto2(t) + out{n,i,s}.rres{t};
    end
end

aerr = toto1 / nbExp;
err = toto2 / nbExp;
figure; plot(aerr);
figure; plot(err(10:end));

% %% Compare the CGPT
% 
% figure; plot(squeeze(Dico.CGPTt{n,s}(1,1,:))/Dico.Scl(s)); hold on; plot(squeeze(CGPTt{n,s}(1,1,:)), 'r')
% toto = Dico.CGPTt{n,s}(1,1,:) - CGPTt{n,s}(1,1,:);
% figure; plot(squeeze(toto(1,1,:)));
% 
% norm(squeeze(toto(1,1,:)), 'fro')/norm(squeeze(Dico.CGPTt{n,s}(1,1,:)),'fro')
% 
