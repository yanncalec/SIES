%% Plot the acquisition systems of limited view in group and the SVD of the operators

%% Add path
clear all;
% close all;
clc;
addpath('~/OOP/');

%% Definition of the small inclusion

delta = 1 ; % diameter of the standard shape
%%
% Initialize an object of |C2boundary|

% B = shape.Ellipse(delta,delta/2,[0,0]',0,2^10);

B = shape.Flower(delta/2,delta/2,[0,0]',0,2^10,5,0.4,0); 
% B = shape.Triangle(delta/2, pi*0.8, 2^10);
% B = shape.Rectangle(delta,0.5*delta,[0,0]',0,2^10);
% % or load image from file
% B = shape.Imgshape('../images/Letters/R.png', delta, delta, 2^10);
% figure; plot(B); axis image;

%%
% Set the conductivity and the permittivity of the inclusion
B.cnd = 3; B.pmtt = 3; B.pmeb = 3;

%% Set up an environment for experience
% The sources/receptors are distributed on a circle whose center is closed
% to the mass of center of the inclusion, up to a small offset.

%%
% center of the measurement circle
mcenter = [0,0]';
%%
% radius of the measurement circle
mradius = 3;

%%
% Make the acquisition configuration. With the class |acq.Planewave|, we make equally distributed
% sources/receivers on a circle of center |mcenter| with radius |mradius|.
N0 = 20; % Number of sources (and receivers)
cfg = acq.Group_Planewave(N0, mcenter, mradius, N0, [5, 0.3*pi, 2*pi]);

%%
% Scaling, translation and rotation parameters:
% scl = 1.5; trl = [-0.25,0.25]'; rtn = pi/3;
scl = 1; trl = 2*[-0.25,0.25]'; rtn = pi/3;
% scl = 1; trl = [0,0]'; rtn = 0;
% trl = D.diameter/2*(rand(2,1)-0.5); 

%%
% The true inclusion _D_ is the shape _B_ after some rigid transform and
% scaling. Apply these transforms:
D = (B < rtn) * scl + trl; % or D = B.t_s_r(trl, scl, rtn);

%%
freq = 2*pi; % Working frequency of sources
pmtt_bg = 1; 
pmeb_bg = 1;

P = PDE.Helmholtz_R2(D, cfg, freq, pmtt_bg, pmeb_bg); 
fig = figure; plot(P); axis image; 
saveas(fig, ['../figures/Acqsys_limview_group.eps'], 'psc2');

%% SVD of operator

Lm = {}; S={};

for n=1:3
    ord = n*5;
    Op = PDE.Helmholtz_R2.make_linop_SCT(cfg, P.wavenb_bg, ord); % Construct the linear operator L

    idx = 1;
    Lm{n} = zeros(cfg.Ns*cfg.Nr, (2*ord+1)^2);

    for rr=1:2*ord+1
        for cc=1:2*ord+1
            CB=zeros(2*ord+1);
            CB(rr,cc)=1;
            toto = Op.L(CB(:), 'notransp');
            Lm{n}(:, idx) = toto;
            idx = idx+1;
        end
    end

    S{n} = svd(Lm{n});
end

fig1=figure; 
for n=1:3
    semilogy(sort(S{n}, 'descend'));
    hold on; grid on;
end

% figure;
% semilogy(sort(S{1}, 'descend'));     hold on; grid on;
% semilogy(sort(S{2}, 'descend'), '--');
% semilogy(sort(S{3}, 'descend'), 'r-.');

saveas(fig1, '../figures/SVD_limview_group.eps','psc2');
