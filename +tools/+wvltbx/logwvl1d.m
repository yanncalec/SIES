clear all
close all

%% Compute the inner product <log, psi> by quadrature
wname = 'db4';

[Tphi, Tpsi, Tval] = wavefun(wname, 12);
[N1, N2] = wavsupport(wname); wvlsize=N2-N1;
Nw = length(Tval); dh = wvlsize/Nw;

j=-1; 

rx=10; ry=11;
nmin0=floor(max(1, -N1+1));
nmin = max(nmin0, ceil(2^(-j)*rx-N2)) 
nmax = floor(2^(-j)*ry-N1)

xs = 2^j*linspace(rx, ry, Nw);
% figure; plot(xs, log(xs)); 

acoeff=[]; dcoeff=[]; 
tt = linspace(N1, N2, Nw);

for n=nmin:nmax
    Y=log(2^j*(tt+n))*2^(j/2);
    acoeff(n-nmin+1) = sum(Tphi.*Y)*dh;
    dcoeff(n-nmin+1) = sum(Tpsi.*Y)*dh;
end
figure; plot(acoeff); 
figure; plot(dcoeff); 

%% Compute the coefficient on a ROI by FWT
N0=1024*16;
sx = linspace(rx, ry, N0);
N1=N0/16;
nlvl = wmaxlev(floor((N0-N1)/2), wname)
idx=wkeep(1:N0, N1);
% nlvl = 8;

[mask, maskc] = wmask_ROI(wname, N0, N1, nlvl);

X=log(sx);
[C, L0] = wavedec(X, nlvl, wname);
W = wvl_vec2scl(C, L0);
Wa = extract_active_coeffs(W, maskc);

for n = 1:length(W)    
    figure; 
    subplot 121; plot(Wa{n});
    subplot 122; plot(W{n}); 
end

Walow = extract_lowscls(Wa,0);
Ca = cell2mat(Walow);

Xr = waverec(Ca, L0, wname);
Err = Xr-X; 
Err = wkeep(Err, N1); % Err=Err(idx);
norm(Err)
figure; plot(Err);
figure; plot(Xr); hold on; plot(X, 'r');
figure; plot(Xr-X);

