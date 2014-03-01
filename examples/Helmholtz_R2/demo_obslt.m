%% Verify the forward operator and translation formula of SCT

%% Add path
clear all;
close all;
clc;
addpath('~/OOP/');

%% Definition of the small inclusion

delta = 1 ; % diameter of the standard shape
%%
% Initialize an object of |C2boundary|

% B = shape.Ellipse(delta,delta/2,[0,0]',0,2^10); % 
B = shape.Flower(delta/2,delta/2,[0,0]',0,2^10,5,0.4,0); 
% B = shape.Triangle(delta/2, pi*0.8, 2^10);
% B = shape.Rectangle(delta,0.5*delta,[0,0]',0,2^10);
% % or load image from file
% B = shape.Imgshape('../images/Letters/R.png', delta, delta, 2^10);
% figure; plot(B); axis image;

%%
% Set the conductivity and the permittivity of the inclusion
B.cnd = 3; B.pmtt = 3; B.pmeb = 3;

%%
% Scaling, translation and rotation parameters:
% scl = .5; trl = [-0.5,0.5]'; rtn = pi/3;
scl = 1; trl = [0,0]'; rtn = pi/3;

%%
% The true inclusion _D_ is the shape _B_ after some rigid transform and
% scaling. Apply these transforms:
D = (B < rtn) * scl + trl; % or D = B.t_s_r(trl, scl, rtn);
figure; plot(D); axis image;

%% Set up an environment for experience
% The sources/receptors are distributed on a circle whose center is closed
% to the mass of center of the inclusion, up to a small offset.

%%
% offset of the measurement center
% vh = D.diameter/2*(rand(2,1)-0.5); 
vh= [0,0]';
%%
% center of the measurement circle
mcenter = D.center_of_mass + vh; 
%%
% radius of the measurement circle
mradius = D.diameter*2; 

%%
% Make the acquisition configuration. With the class |acq.Planewave|, we make equally distributed
% sources/receivers on a circle of center |mcenter| with radius |mradius|.

cfg = acq.Planewave(5, mcenter, mradius, 10, [1, 2*pi, 2*pi]);
% figure;plot(cfg); axis image;

%%
freq = 5; % Working frequency of sources
pmtt_bg = 1; 
pmeb_bg = 1;
k0 = sqrt(pmtt_bg*pmeb_bg)*freq;

P = PDE.Helmholtz_R2(D, cfg, freq, pmtt_bg, pmeb_bg); 

figure; plot(P, 'LineWidth', 1); axis image;

% xs = cfg.src(1); Sx = linspace(-3,3,100); 
% F = PDE.Conductivity_R2.solve_forward(D, 0, xs, Sx, Sx);
% figure; imagesc(F); colorbar()

%% Simulation of the MSR data
% freqlist = linspace(0,1,10);
tic
data = P.data_simulation();
toc

%%
% % Calculate the field and plot it. 
% [F, F_bg, SX, SY, mask] = P.calculate_field(3, cfg.src(1), 2*D.diameter, 64);
% P.plot_field(F, F_bg, SX, SY, 10, 'LineWidth', 2);

%%
% add white noise
nlvl = 0.1; 
data = P.add_white_noise(data, nlvl);

%%
% Don't forget to take the transpose, each row needs to correspond to a source
MSR = data.MSR{1}.'; 
% MSR = data.MSR_noisy.'; 

%% Compute the SCT
%%
% maximum order of the reconstruction
ord = 10; 

%%
% The inclusion corresponding to the reconstruction
E = D - mcenter;

tic
W0 = D.SCT(ord, freq, pmtt_bg, pmeb_bg); 
toc
figure; imagesc(abs(W0)); colorbar()

%% test the forward operator
tic
out = P.make_linop_SCT(cfg, k0, 10);
toc
Op = out.L;

Y = Op(W0, 'notransp'); 
norm(Y-MSR(:))/norm(MSR(:))

Ym = reshape(Y, cfg.Ns, cfg.Nr); Ym = Ym.';
Y0=MSR.';

figure; 
subplot(221); plot(real(Ym(:))); title('Forward operator: real part');
subplot(222); plot(imag(Ym(:))); title('Forward operator: imag part');
subplot(223); plot(real(Y0(:))); title('Simulated data: real part');
subplot(224); plot(imag(Y0(:))); title('Simulated data: imag pat');

% figure; 
% plot(real(Ym(:))); hold on; plot(real(Y0(:)), 'r'); title('real part');
% figure; 
% plot(imag(Ym(:))); hold on; plot(imag(Y0(:)), 'r'); title('imag part');

err = Ym-Y0;
figure; 
subplot(121); plot(real(err(:))); title('Error: real part');
subplot(122); plot(imag(err(:))); title('Error: imag part');

%% Verify the translation formula
W0=D.SCT(50, 1, 1, 1);

z0 = [2,2]';
E = D+z0;
W1=E.SCT(8, 1, 1, 1);

L=10; G=zeros(2*L+1,1);
for n=-L:L
    G(n+L+1)=tools.Helmholtz.cylind_wave(1, n, z0);
end

% W2c = conv2(conj(G), G, W0); 
W2c = conv2(G, conj(G), W0); 
W2=wkeep(W2c, size(W1));
figure; imagesc(abs(W1-W2)); colorbar();
norm(W1-W2,'fro')

%%
% With the translation formula in Lemma 5.3 the result is wrong:

W3=zeros(size(W1));
for r=1:size(W1,1)
    for c=1:size(W1,2)
        m = r-1-(size(W1,1)-1)/2;
        n = c-1-(size(W1,2)-1)/2;
        for l=-L:L
            
            if m+l>=-(size(W0,1)-1)/2 && m+l<=(size(W0,1)-1)/2
                if n+l>=-(size(W0,2)-1)/2 && n+l<=(size(W0,2)-1)/2                
                    rr = m+l+(size(W0,1)-1)/2+1;
                    cc = n+l+(size(W0,2)-1)/2+1;
                    W3(r,c) = W3(r,c) + abs(G(l+L+1))^2 * W0(rr, cc);
                end
            end
        end
    end
end

figure; imagesc(abs(W3)); colorbar();
figure; imagesc(abs(W1-W3)); colorbar();
norm(W1-W3,'fro')
