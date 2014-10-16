function A = make_matrix_WPT(Xs, Z, WF, nbScl, stdmode)
% A = make_matrix_WPT(Xs, Z, WF, nbScl, stdmode)
% Construct the matrix (the linear operator L) involved in the model
% Conductivity_R2 with wavelet
% Inputs:
% Xs: coordinates of sources/receivers, of dimension 2X?
% Z: reference center
% WF: an object of wavelet.Frame class
% nbScl: number of scales invovled in the computation of WPT
% stdmode: true for standard representation

if nargin < 5
	stdmode = 1;
end
if nargin < 4
	nbScl = 0;
end

% Xs = reshape(Xs, 2, []);
N = size(Xs,2); % number of sources

mask0 = WF.mask_active_lowscls(nbScl); % mask for active coefficients
mask = WF.invert_scl_v(mask0, nbScl);
% [Am, Dm] = WF.vec2scl(mask0);
% size(Dm)
% mask = invert_scl(Am, Dm(1:nbScl,:)); % invert the scale order
aidx = find(mask); % index of active wavelets

if stdmode
	% A = zeros(N, WF.nbWvl);
	% A = spalloc(N, WF.nbWvl, WF.nbApprx*N);
	
	A = zeros(N, length(aidx));
	
	for s=1:N
		if mod(s-1,10)==0
			fprintf('Processing the source/receiver No. %d\n', s);
		end
		
		C = WF.decomp_Green2D(Xs(:,s)-Z);
		% C = WF.decomp_Green2D(Xs(:,s));
		C = WF.invert_scl_v(C, nbScl);
		% [Ac, Dc] = WF.vec2scl(C);
		% C = invert_scl(Ac, Dc(1:nbScl,:));
		A(s,:) = reshape(C(aidx), 1, []);
	end
else
	% A = Conductivity_R2.make_matrix_WPT_nonstd(Xs, Z0, WF);
	error('Not implemented');
end
end
