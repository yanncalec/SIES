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

N0 = 40; % Number of sources (and receivers)
pmtt_bg = 1; pmeb_bg = 1;

cfg = acq.Group_Planewave(N0, [0,0]', 3, N0, [4, 1*pi, 2*pi]);
P = PDE.Helmholtz_R2(D, cfg, freq, pmtt_bg, pmeb_bg); 
figure; plot(P, 'LineWidth', 1); axis image;

data = P.data_simulation();
MSR = data.MSR{1};

ord0 = 30;
W0 = D.SCT(ord0, freq, pmtt_bg, pmeb_bg);

%% SCT and fwd operator

% ord = ord0; W = W0;
% Op = PDE.Helmholtz_R2.make_linop_SCT(P.cfg, P.wavenb_bg, ord);

% derr = Op.L(W,'notransp') - MSR(:);
% norm(derr) / norm(MSR,'fro')

%% recon
% ord = 20;
% W = wkeep(W0, [2*ord+1, 2*ord+1]);
% out = P.reconstruct_SCT(data.MSR{1}, ord, 1e5, 1e-4);    
% err = norm(out.SCT - W , 'fro') / norm(W, 'fro')

%% with noise
ordmin = 8;
ordmax = 15; %ordmin+floor((min(P.cfg.Ns, P.cfg.Nr)-1)/2);

Nlvl = 3;
Nexp = 50;

err = zeros(ordmax-ordmin+1, Nlvl, Nexp);
rerr = err;

for ord=ordmin:ordmax
    W = wkeep(W0, [2*ord+1, 2*ord+1]);
    oidx = ord-ordmin+1;
    
    for n=1:Nlvl
        nlvl = n/(Nlvl+1) * 0.1

        for m=1:Nexp
            data = P.add_white_noise(data, nlvl);
            toto = P.reconstruct_SCT(data.MSR_noisy{1}, ord, 1e5, 1e-4);    
            err(oidx,n,m) = norm(toto.SCT - W , 'fro');
            rerr(oidx,n,m) = err(oidx,n,m) / norm(W, 'fro');
        end
    end
end

Rerr = squeeze(mean(rerr,3));

% E = Rerr(1:end, :); 

fig1 = figure; hold on; grid on;
for n=1:Nlvl-1
    plot(ordmin:ordmax, Rerr(:,n), 'LineWidth', 1); 
end

saveas(fig1, '../figures/res_ord_limview.eps', 'psc2');
