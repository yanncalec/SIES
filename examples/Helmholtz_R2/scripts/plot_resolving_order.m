%% Add path
clear all;
close all;
clc;
addpath('~/OOP/');

nbPoints = 2^10;
delta = 1;

B = shape.Flower(delta/2,delta/2,[0,0]',0,nbPoints,5,0.4,0);
B.cnd = 3;    
B.pmtt = 3;
B.pmeb = 3;

%% data simulation
freq = 2*pi;

scl = 1; trl = [0,0]'; rtn = pi/3;
D = (B<rtn) * scl + trl;

N0 = 91;
pmtt_bg = 1; pmeb_bg = 1;

cfg = acq.Planewave(N0, [0,0]', 3, N0, [1, 2*pi, 2*pi]);
P = PDE.Helmholtz_R2(D, cfg, freq, pmtt_bg, pmeb_bg); 
figure; plot(P, 'LineWidth', 1); axis image;

data = P.data_simulation();
MSR = data.MSR{1};

ord0 = floor((min(P.cfg.Ns, P.cfg.Nr)-1)/2);
W0 = D.SCT(ord0, freq, pmtt_bg, pmeb_bg);

%% SCT and fwd operator

% ord = 30;
% W = wkeep(W0, [2*ord+1,2*ord+1]);

% Op = PDE.Helmholtz_R2.make_linop_SCT(P.cfg, P.wavenb_bg, ord);

% derr = Op.L(W,'notransp') - MSR(:);
% norm(derr) / norm(MSR,'fro')

%% recon

ord = floor((min(P.cfg.Ns, P.cfg.Nr)-1)/2);
W = wkeep(W0, [2*ord+1,2*ord+1]);

out = {};
Nlvl = 5;
Nexp = 100;

err = zeros(ord+1, Nlvl, Nexp);

for n=1:Nlvl
    nlvl = n/Nlvl

    for m=1:Nexp
        data = P.add_white_noise(data, nlvl);
        toto = P.reconstruct_SCT_analytic(data.MSR_noisy{1}, ord);    
        bobo = tools.matrix_norm_by_order(toto.SCT - W, 'center');
        err(:,n,m) = bobo(:);
    end
end

Err = squeeze(mean(err,3));
uouo = tools.matrix_norm_by_order(W, 'center');
Rerr = Err ./ repmat(uouo(:), 1, Nlvl);

fig1 = figure; hold on;
for n=1:Nlvl
    plot(Rerr(:,n), 'LineWidth', 1); 
end
xlim([0,50]); ylim([0,0.7]); grid on;
plot(0.1*ones(ord+1,1), 'r', 'LineWidth',1);

saveas(fig1, '../figures/res_ord_fullview.eps', 'psc2');
% plot(Rerr(:,2), '-.');
% plot(Rerr(:,3), '--');
% plot(Rerr(:,4), 'r-.');
% plot(Rerr(:,5), 'r-.');

