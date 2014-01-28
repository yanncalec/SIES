%% Demo of the wavelet.OrthoWvl class
% This script shows how to use |PDE.Conductivity_R2| class for data simulation
% and |wavelet.OrthoWvl| class for reconstruction of WPTs.

%% Add path
clear all;
close all;
clc;
addpath('../../');
% matlabpool;

%% Definition of the small inclusion

delta = 1 ; % diameter of the standard shape
%%
% Initialize an object of |C2boundary|

B = shape.Flower(delta/2,delta/2,[0,0]',0,2^12,5,0.4,0); 
% figure; plot(B); axis image;
cc
%%
% Set the conductivity and the permittivity of the inclusion
B.cnd = 4/3; B.pmtt = 1;

%%
% Scaling, translation and rotation parameters:
% scl = .5; trl = [-0.5,0.5]'; rtn = pi/3;
scl = 1; trl = [0,0]'; rtn = 0;

%%
% The true inclusion _D_ is the shape _B_ after some rigid transform and
% scaling. Apply these transforms:
D = (B < rtn) * scl + trl; % or D = B.t_s_r(trl, scl, rtn);
% figure; plot(D); axis image;

%% Construction of Daubechies wavelet basis
% clear WF
ROI_width = D.diameter*1.1;
Rsl = ROI_width/16;
disp('Constructing the wavelet object...');
tic
WF=wavelet.OrthoWvl('Daubechies', 6, D.center_of_mass, ROI_width, Rsl);
toc

%%
% Theoretical values of WPT
nbScl = 0;
disp('Computing theoretical WPT...');
tic
WPT = D.WPT_std(WF, 0, nbScl);
toc

%%
% Imaging
% WPT11 = WPT.onedir(1,1);
% figure; imagesc(abs(WPT11.data)); axis tight; colorbar()
% figure; imagesc(log10(abs(WPT11.data))); axis tight; colorbar()

% figure; semilogy(sort(abs(WPT.AA(:)),'descend')); axis tight;

% DD = WPT.DD{3,3}{1,1};
% Wim0 = WPT.diag_imaging(DD);
% figure; imagesc(abs(Wim0)); axis tight; colorbar; title('Imaging of X');
% figure; imagesc(log10(abs(Wim0))); axis tight; colorbar; title('Imaging of X in log');
% figure; semilogy(sort(abs(DD(:)),'descend')); axis tight;

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
E = D-mcenter;
%%
% radius of the measurement circle
mradius = ROI_width/sqrt(2) * 1.1; 
% mradius = WF.AROI.width/sqrt(2)*1.1;

%%
% Make the acquisition configuration. With the class
% |acq.Coincided_circle|, we make equally distributed sources/receivers on
% a circle of center |mcenter| with radius |mradius|.

cfg = acq.Coincided_circle(mcenter, mradius, 200, [1, 2*pi, 2*pi]);
% cfg = acq.Coincided_circle(mcenter, mradius, 10, [5, 0.2*pi, 2*pi]);
% cfg = acq.Coincided_circle(mcenter, mradius, 50, [1, pi, 1*pi]);

%%
% With the class |acq.Group_circle|, the sources/receivers are divided into
% groups and they communicate only within the same group.

% cfg = acq.Group_circle(mcenter, mradius, 20, 21, 0.3*pi, 2*pi); 
% cfg = acq.Group_circle(mcenter, mradius, 10, 5, 0.2*pi, 1*pi);

%%
% clear P;
freq = 0; % Working frequency of sources
P = PDE.Conductivity_R2(D, cfg, freq); 

figure; plot(P, 'LineWidth', 1); axis image;

% xs = cfg.src(1); Sx = linspace(-3,3,100); 
% F = PDE.Conductivity_R2.solve_forward(D, 0, xs, Sx, Sx);
% figure; imagesc(F); colorbar()

%% Simulation of the MSR data
% freqlist = linspace(0,1,10);

disp('Simulation of data...');
tic
data = P.data_simulation();
toc

%%
% add white noise
nlvl = 0.1; 
data = P.add_white_noise(data, nlvl);

%%
% Don't forget to take the transpose, each row needs to correspond to a source
MSR = data.MSR.'; 
% MSR = data.MSR_noisy.'; 

% %%
% % Calculate the field and plot it. 
% xlim = mcenter(1) + [-1,1]*D.diameter*3; ylim = mcenter(2) + [-1,1]*D.diameter*3; dh = 0.05;
% [F, F_bg, SX, SY] = P.calculate_field(1, xlim, ylim, dh);
% P.plot_field(1, F, F_bg, SX, SY, 25);

%% Make the forward operator and verify the error
disp('Making linear operators...');
tic
Op = P.make_linop_WPT(cfg, WF, nbScl);
toc

Y = Op.L(WPT.data, 'notransp');
err = MSR(:) - Y(:);
fprintf('L2 error and relative error of forward operator: %e, %e\n', norm(err), norm(err)/norm(MSR(:)));
fprintf('Linf error of forward operator: %e\n', max(abs(err)));
figure; plot(Y(:)); hold on; plot(MSR(:), 'r'); title('Data and L(X)');
figure; plot(Y(:)-MSR(:)); title('Error of L(X)');

% Yt = reshape(Y, size(MSR)); Yt=Yt.'; MSRt = MSR.';
% figure; plot(Yt(:)); hold on; plot(MSRt(:), 'r');
% figure; plot(Yt(:)-MSRt(:));

%%
% Verify that L(W)-L(W^t) ~ 0
Yt = Op.L(WPT.data.', 'notransp');
% figure; plot(Yt(:)-MSR(:));
fprintf('L2 error of L(W-W^t): %e\n', norm(Y-Yt));
% figure; plot(Y-Yt);

% but W<>W^t
fprintf('L2 error of W-W^t: %e\n', norm(WPT.data-WPT.data.','fro'));

%% Increasing the number of scales can reduce the error of the fwd operator
% disp('Making linear operators...');
% tic
% Op1 = P.make_linop_WPT(cfg, WF, 1);
% toc

% WPT1 = WPT.lowscls(1);
% Y = Op1.L(WPT1.data, 'notransp');
% err = MSR(:) - Y(:);
% fprintf('L2 error and relative error of forward operator: %e, %e\n', norm(err), norm(err)/norm(MSR(:)));
% fprintf('Linf error of forward operator: %e\n', max(abs(err)));
% figure; plot(Y(:)); hold on; plot(MSR(:), 'r'); title('Data and L(X)');
% figure; plot(Y(:)-MSR(:)); title('Error of L(X)');

%% Approximation error of the acq-sys operator for CGPT
disp('*********CGPT Case*********');

ord = WF.wvlord;
CGPT = E.CGPT(ord, freq);
Opc = P.make_linop_CGPT(cfg, ord, 0);
Y = Opc.L(CGPT.data,'notransp');

err = MSR(:) - Y(:);
fprintf('L2 error and relative error of forward operator: %e, %e\n', norm(err), norm(err)/norm(MSR(:)));
fprintf('Linf error of forward operator: %e\n', max(abs(err)));

figure; plot(Y(:)); hold on; plot(MSR(:), 'r'); title('CGPT: Data and L(X)');
figure; plot(Y(:)-MSR(:)); title('CGPT: Error of L(X)');

%% Calculate the mask
% AA = full(WPT.AA.*mask); tools.l0norm(AA)/numel(AA)

AA =WPT.AA;
Nz = ceil(0.005*numel(AA)); % floor(numel(AA)*0.001); 

[AA1,Err1,~] = tools.Best_Napprx(AA, Nz); 
Y1 = Op.L(AA1, 'notransp');
% Err1/norm(WPT.AA,'fro')
% norm(MSR(:) - Y1(:),'fro')/norm(MSR(:),'fro')

figure; plot(Y1(:)); hold on; plot(MSR(:), 'r'); 
title(['Data and L(Xc) - compressed WPT']);
mask1 = full(double(abs(AA1)>0));

mask2 = WPT.make_WPT_mask(WF.ApprxSpace{1}.dim,3);

sum(sum(abs(mask1.*mask2-mask1)))/sum(sum(mask1))
sum(mask2(:)) / numel(mask2)

figure; imagesc(mask1); axis image; title('Mask for WPT');
figure; imagesc(mask2); axis image; title('Mask for WPT');

Wim = WPT.diag_imaging(mask1);
figure; imagesc(abs(Wim)); axis tight; colorbar; title('Imaging of mask');

%% Reconstruction of WPT
disp('Reconstruct wavelet tensors...');

% L1 reconstruction
tic
out1 = P.reconstruct_WPT_l1(MSR(:), Op.L, mask, 1e-6, 100000, 1e-7, 'fista');
toc

% Least-square reconstruction
% tic
% out2 = P.reconstruct_WPT(MSR(:), Op.L, mask2, 10000, 1e-7);
% toc

Xr=out2.X(:);  
sum(abs(Xr(:))>0)/numel(Xr)
figure; semilogy(sort(abs(Xr), 'descend')); axis tight; hold on; 
semilogy(sort(abs(WPT.AA(:)), 'descend'), 'r'); axis tight;
title('X vs Xr: decay');

Ar = reshape(Xr, WF.nbApprx, WF.nbApprx);

Yr = Op.L(Xr, 'notransp'); 
errr = MSR(:) - Yr(:);
fprintf('L2 error and relative error of reconstruction: %e, %e\n', norm(errr), norm(errr)/norm(MSR(:)));
fprintf('Linf error of reconstruction: %e\n', max(abs(errr)));

figure; plot(Yr(:)); hold on; plot(MSR(:), 'r'); title('Data and L(Xr)');
% figure; plot(Yr(:)-MSR(:));

figure; imagesc(abs(Ar)); axis tight; colorbar; title('Xr');
figure; imagesc(log10(abs(Ar))); axis tight; colorbar; title('log Xr');
% figure; imagesc(abs(WPT.AA)); axis tight; colorbar; title('X');
% figure; imagesc(log10(abs(WPT.AA))); axis tight; colorbar; title('log X');

% Imaging
Wim = WPT.max_imaging(Ar);
% Wim = WPT.diag_imaging(Ar);
figure; imagesc(abs(Wim)); axis tight; colorbar; title('Imaging of Xr');
figure; imagesc(log10(abs(Wim))); axis tight; colorbar; title('Imaging of Xr in log');

% Wim0 = WPT.diag_imaging(WPT.AA);
% figure; imagesc(abs(Wim0)); axis tight; colorbar; title('Imaging of X');
% figure; imagesc(log10(abs(Wim0))); axis tight; colorbar; title('Imaging of X in log');

cc
%% Theoretical values of WPT

% AA = WPT.AA;
% res = norm(MSR(:) - Op.L(AA(:), 'notransp'))/norm(MSR(:)) % very large error?  WF.Rsl^2
% norm(MSR(:))

% %%
% % Imaging by WPT
% Nz = 5000 % floor(numel(AA)*0.001); sqrt(Nz)
% [AA1,Err1,~] = tools.Best_Napprx(AA, Nz);
% Err1/norm(AA,'fro')
% figure; semilogy(sort(abs(AA(:)), 'descend')); axis tight;
% figure; semilogy(sort(abs(AA1(:)), 'descend')); axis tight;

% Aim = WF.WPT_imaging(AA);
% figure;imagesc(log10(abs(Aim)));colorbar(); axis image;
% Aim1 = WF.WPT_imaging(AA1);
% figure;imagesc(log10(abs(Aim1)));colorbar(); axis image;

% NN=1896; % index of the wavelet
% DD_img = reshape(AA(NN,:),WF.ROI_cdim(1,:)); 

% figure; imagesc(abs(DD_img)); colorbar();  axis image;
% figure; imagesc(log10(abs(DD_img))); colorbar(); axis image;
% img=zeros(WF.ROI_cdim(1,:)); img(NN)=1; 
% figure;imagesc(img);colorbar(); axis image;

%% X=pinv(out.As)*MSR*pinv(out.Ar');
Xm = reshape(out.X, WF.nbApprx, WF.nbApprx);
[uu,vv]= max(abs(Xm));
figure;plot(uu)

G = zeros(WF.nbApprx,1);
for n=1:WF.nbApprx    
    [uu, vv] = max(abs(Xm(n,:)));
    G(vv) = G(vv)+Xm(n,vv);
end
figure; imagesc(reshape(G, WF.ROI_cdim(1,:))); colorbar

M1 = reshape(Xm(:, 10), WF.ROI_cdim(1,:));
figure; imagesc(log10(abs(M1))); colorbar

%% 
%
DD=cell2mat(out.WPT.DD{1,1});
figure;
imagesc(log10(abs(DD))); colorbar()

% Verify the sparsity of wavelet tensors by keeping 5% of the largest non
% zeros coefficients
N = ceil(numel(DD) * 0.05); 
[DD1, err, ~] = tools.Best_Napprx(DD, N);
% Relative error of approximation
err/norm(DD,'fro')
figure; imagesc(log10(abs(DD1))); colorbar(); 

NN=110; % index of the wavelet
toto=out.WPT.DD{1,1}{1,1}; % Direction 1-1, scale 1-1
DD_img = reshape(toto(NN,:),WF.dmesh.shape{1}); 
figure; imagesc((abs(DD_img))); colorbar(); 
figure; imagesc(log10(abs(DD_img))); colorbar(); 
% img=zeros(WF.dmesh.shape{1}); img(NN)=1; 
% figure;imagesc(img);colorbar();

%% Study of the (lambda*I-Kstar)^-1 operator
Kstar=ops.Kstar(D,'P0',1);
Id=ops.Ident(D,'P0',1);
lambda=0.5+1e-5;
A=lambda*Id.stiffmat-Kstar.stiffmat;
% A=0.5*Id.Kmat-Kstar.Kmat;
Ai=inv(A); cond(A)
% figure; imagesc(abs(Kstar.stiffmat)); axis tight; colorbar
figure; imagesc(log10(abs(Ai))); axis tight; colorbar
figure; imagesc(log10(abs(A))); axis tight; colorbar

N = 100; figure; plot(Ai(:,N)); 

cc
%% plot a wavelet
% N=256;
% ss=linspace(-5,5,N);
% [SX, SY] = meshgrid(-ss, ss); X=[SX(:) SY(:)]';
% [val, dX, dY] = wavelet.Gabor.eval(X, [0,0]', 0, 1, parms.pwLap, parms.omega);
% % [val, dX, dY] = wavelet.Gabor.eval(X, [0,0]', 0, 1, 0, [1e0, 0]);
% V=real(reshape(val,N,N)); %V=V.*(V>1e-2);
% figure; surf(SX,SY,V); colorbar()
% figure; imagesc(V); colorbar()

% val1 = exp(-ss.^2/2 + 2*pi*1i*ss*10);
% figure;plot(real(val1));

%%

% AA=reshape(W.AA(18,:),apprx.shape);
% figure;imagesc((abs(AA))); colorbar(); 
% 
% figure; plot(B); axis equal; hold on; plot(apprx.pos(1,:),apprx.pos(2,:),'.');
% 
% figure; plot(B); axis equal; hold on; plot(detail.pos{1}(1,:),detail.pos{1}(2,:),'.');
% figure; plot(B); axis equal; hold on; plot(detail.pos{2}(1,:),detail.pos{2}(2,:),'.');
% figure; plot(B); axis equal; hold on; plot(detail.pos{3}(1,:),detail.pos{3}(2,:),'.');

% S = zeros(prod(WF.dmesh.shape{sidx}),1);
% T=zeros(size(toto,1),1);
% 
% for nn=1:size(toto,1)
%     [~,idx]=max(abs(toto(nn,:)));
%     T(nn) = idx;
%     S(idx) = S(idx) + 1;
% end



% %% a
% cc
% figure;imagesc(log10(abs(DD_img))); colorbar(); axis image;
% img=zeros(WF.dmesh.shape{1}); img(NN)=1; 
% figure;imagesc(img); axis image; colorbar();

% % [T, ~] = max(abs(toto), [], 2); 
% % S_img = reshape(T,WF.dmesh.shape{sidx});
% % figure; imagesc(S_img.*(S_img>10)); colorbar(); axis image;
% % figure; imagesc(log10(S_img)); colorbar(); axis image;

% % %figure;imagesc(log10(abs(DD_img))); colorbar(); axis image;

% L=Op.L;
% save('linsys.mat', 'L'); , 'cfg', 'P', 'MSR');

% for n=1:numel(WPT.AA)    
%     zz=WPT.AA(n,:);
%     [~, idx(n)]=max(abs(zz));
% end

% mm=525;
% zz=WPT.AA(:,mm);
% [~, idx]=max(abs(zz)); idx

% [xx,yy]=sort(abs(zz), 'descend');
% xx=full(xx); nn=100;
% uu=zeros(36^2,1);
% uu(yy(1:nn)) = zz(yy(1:nn));
% figure;imagesc(reshape(uu,36,36)); axis tight; colorbar

%% test
% A=rand(3,5); B=rand(4,5);
Q = @(X) reshape(A * reshape(X, 5, 5) * B', [] ,1);
X=rand(5,5);
Q(X)
Qm = zeros(12, 25);
for m=1:12
    for n=1:25
        E = zeros(25,1); 
        E(n) = 1;
        Qm(:,n) = Q(E);
    end
end
Qmn = tools.row_normalization(Qm);
An = tools.row_normalization(A);
Bn = tools.row_normalization(B);

Qn = @(X) reshape(An * reshape(X, 5, 5) * Bn', [] ,1);
Qn(X) - Qmn*X(:)

%% L1-Fista test
N = 100;
A = randn(N); 
X = randn(N,1); X = tools.Best_Napprx(X, 15);
Y = A*X;
out = tools.solver.l1fista(A, Y, [], 1e1, 1e-10, 1000, 1);

figure;plot(out.X); hold on; plot(X,'r')



% Sr = sort(abs(WPT.AA(:)), 'descend'); Sr=Sr(find(Sr>0));
% figure; semilogy(Sr); axis tight;

% yy=(1:length(Sr)).^(-1/0.9838);
% yy=10.^(-2*10^-5*(1:length(Sr)));
% figure;semilogy(yy); axis tight
