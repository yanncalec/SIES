%% Matching in a dictionary

clear all;
addpath('~/SIES/');

%% Load the dictionary

%pathname = '~/Data/measurements/Pulse/Transformed/0.125pi/';
% pathname = '/Volumes/Yue/Pulse/Original/2pi/';
pathname = '/Volumes/Yue/Pulse2/Original/2pi/';
load([pathname, 'data8_6scl.mat']);

%% Build shape descriptors
load ~/Data/dico/Pulse/smalldico3_21scl.mat;

SD_method = 2;

for m=1:length(Dico.B) % iteration on the shape
    for s = 1:scl        
        [Dico.CGPTt{m,s}, Dico.dt(s)] = asymp.CGPT.CGPT_time_truncation(Dico.CGPTt0{m,s}, Dico.dt0(s), Data.Tmax(s), Data.Ntime);
    end

%     toto = dico.CGPT.ShapeDescriptor_PT_time(Dico.CGPTt(m,:), Dico.Scl, SD_method);
%     Dico.SDt{m} = toto(Dico.extrema, :);

    Dico.SDt{m} = dico.CGPT.ShapeDescriptor_PT_time(Dico.CGPTt(m,:), Dico.Scl, SD_method);
end

%% Make sure that the dictionary is good for shape identification
Err = zeros(length(Dico.B), length(Dico.B)); Idx = Err;
for m = 1:length(Dico.B)
    [Err(m,:), Idx(m,:)] = dico.CGPT.SD_Matching_time(Dico.SDt{m}, Dico.SDt, 2);
end

fig1= figure;
bar(Err, 'facecolor', 'none'); 
set(gca, 'XTickLabel', Dico.names, 'XTick',1:length(Dico.B));
hold on; 

% Identification result is plotted in green, if the result is wrong, the
% true shape is marked in red.
bar(eye(length(Dico.B)).*Err, 'r'); 

toto=eye(length(Dico.B)); 
idx = Idx(:,1); 
bar(toto(idx, :).*Err, 'g'); 

%% Plot SD
bidx = [1, 2, 3];
sidx = (-5:-1)+11 ;
SDt = Dico.SDt(bidx);

fig=figure;
wd = 1.;
toto = reshape(SDt{1}(:,sidx), 1, []); plot(toto,'r:', 'LineWidth',wd); hold on;
toto = reshape(SDt{2}(:,sidx), 1, []); plot(toto,'r', 'LineWidth',wd); hold on;
toto = reshape(SDt{3}(:,sidx), 1, []); plot(toto,'g--', 'LineWidth',wd); hold on;
% toto = reshape(SDt{4}(:,sidx), 1, []); plot(toto,'LineWidth',wd); hold on;

legend(Dico.names{bidx});
saveas(fig, '~/features.eps','psc2');

% Compare each shape with the disk
W = SDt{2}(:,sidx); 

fig = figure;
toto = reshape(SDt{2}(:,sidx)-W, 1, []); plot(toto,'r', 'LineWidth',wd); hold on;
toto = reshape(SDt{3}(:,sidx)-W, 1, []); plot(toto,'g--', 'LineWidth',wd); hold on;
% toto = reshape(SDt{4}(:,sidx)-W, 1, []); plot(toto,'LineWidth',wd); hold on;

legend(Dico.names{bidx(2:end)});
saveas(fig, '~/features_comp.eps','psc2');

%% Parameters 
nbShapes = length(Data.B);
Sidx = 1:nbShapes;

Ntime = length(Data.data{1,1}.MSR);

cfg = Data.cfg;
P = PDE.PulseImaging_R2(Dico.B{1}, Dico.cnd, Dico.pmtt, Dico.waveform(2,:), Dico.dt(2), cfg);
figure; plot(P); axis image;

%% Dico-matching
nlvl = 0.5; % noise level
nbExp = 1;

Err = zeros(nbShapes, length(Dico.B)); Idx = Err;
% Err = zeros(nbShapes,nbShapes); Idx = Err;

scl = 6; %length(Dico.Scl);
CGPTt = cell(nbShapes, scl);
out = cell(nbShapes, nbExp, scl);

for m = 1:nbShapes
    errtmp = zeros(nbExp, length(Dico.B));
    %idx = zeros(nbExp, length(B));
    
    for i = 1:nbExp
        for s = 1:scl
            data_noisy = P.add_white_noise(Data.data{m,s}, nlvl);
            % out = P.reconstruct_CGPT(data_noisy.MSR_noisy, 1, 10^5, 10^-5, 1, 'lsqr');
            out{m,i,s} = P.reconstruct_CGPT(data_noisy.MSR_noisy, 1); % pinv
            
            CGPTt{m,s} = tools.cell2mat3D(out{m,i,s}.CGPT);
        end
        
        SDt = dico.CGPT.ShapeDescriptor_PT_time(CGPTt(m,:), Dico.Scl, SD_method);
        % SDt = SDt(Dico.extrema, :);
        [errtmp(i,:), ~] = dico.CGPT.SD_Matching_time(SDt, Dico.SDt, 1:scl);
    end
    
    Err(m,:) = mean(errtmp,1);
    [~, Idx(m,:)] = sort(Err(m,:));
end

% Idx
% Err

% Interpretation of the result
% We show in a bar figure the similarity between dictionary shape descriptors and the one reconstructed from data.

fig1= figure; 
bar(Err, 'facecolor', 'none'); 
set(gca, 'XTickLabel', Dico.names, 'XTick',1:nbShapes);
hold on; 

% Identification result is plotted in green, if the result is wrong, the
% true shape is marked in red.
bar(eye(nbShapes).*Err, 'r'); 

toto=eye(nbShapes); 
idx = Idx(:,1); 
bar(toto(idx, :).*Err, 'g'); 
cc
%% Compare the MSR
n=6;
s=2; 

toto1 = zeros(1, Ntime);
toto2 = zeros(1, Ntime);

for i = 1:nbExp
    for t=1:Ntime
        toto1(t) = toto1(t) + out{n,i,s}.res{t};
        toto2(t) = toto2(t) + out{n,i,s}.rres{t};
    end
end

aerr = toto1 / nbExp;
err = toto2 / nbExp;
figure; plot(aerr);
figure; plot(err(10:end));

%% Compare the CGPT

figure; plot(squeeze(Dico.CGPTt{n,s}(1,1,:))); hold on; plot(squeeze(CGPTt{n,s}(1,1,:)), 'r')
toto = Dico.CGPTt{n,s}(1,1,:) - CGPTt{n,s}(1,1,:);
figure; plot(squeeze(toto(1,1,:)));

norm(squeeze(toto(1,1,:)), 'fro')/norm(squeeze(Dico.CGPTt{n,s}(1,1,:)),'fro')

cc
%% CGPTtF={};

n=3;
for s = 1:scl
    % pulse waveform at the scale s
    % [CGPTtF{s}, ~] = asymp.CGPT.theoretical_CGPT_time_recon((Dico.B{n}<0.2*pi)*2+[1,1]', Dico.cnd, Dico.pmtt, 1, Dico.Tmax0, Dico.Ntime, Dico.Scl(s));
    [CGPTtF{s}, ~] = asymp.CGPT.theoretical_CGPT_time_recon(Dico.B{n}, Dico.cnd, Dico.pmtt, 1, Dico.Tmax0, Dico.Ntime, Dico.Scl(s));
end 
   
% 
% %toto = Dico.CGPTt{3,1};
% toto = Dico.CGPTt{3,1}-CGPTtF{1};

% CGPTtF = CGPTt;

SDt = dico.CGPT.ShapeDescriptor_PT_time(CGPTtF, Dico.Scl);
% SDt = Dico.SDt(n,:);

toto = SDt - Dico.SDt(n,:);
norm(toto)
figure; plot(toto)
figure; plot(Dico.SDt(n,:)); hold on; plot(SDt, 'r')

toto = SDt - Dico.SDt(2,:);
norm(toto)
figure; plot(toto)
figure; plot(Dico.SDt(2,:)); hold on; plot(SDt, 'r')





% %%
% Y0 = tools.cell2mat3D(data_noisy.MSR);
% Y1 = tools.cell2mat3D(data_noisy.MSR_noisy);
% rr=1; cc=1;
% figure; 
% plot(squeeze(Y0(rr,cc,:))); hold on; 
% plot(squeeze(Y1(rr,cc,:)),'r'); 
