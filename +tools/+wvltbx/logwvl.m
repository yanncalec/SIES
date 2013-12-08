clear all;

wname = 'db6';
sx = linspace(3,5, 1024);

[Sx, Sy] = meshgrid(sx, sx);
meshpoints = [Sx(:) Sy(:)]';
X = 1/2*log(Sx.^2+Sy.^2);
reshape(X, size(Sx));
figure; imagesc(X); axis image; colorbar
nlvl = wmaxlev(length(sx), wname)

[C, L] = wavedec2(X, nlvl, wname);
P= L(2:end-1,:); P=repmat(P, 1,3);
P=reshape(P,[], 2);

cc

