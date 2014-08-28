function CGPTt = addnoise_recon(obj, data, nlvl, ord, maxiter, tol, symmode, method, op)
% method = 'pinv' or 'lsqr'

if nargin < 9
    op = [];
end

if nargin < 8
    method = 'lsqr';
end

if nargin < 7
    symmode = 1;
end

if nargin < 6
    tol = 1e-8;
end

if nargin < 5
    maxiter = 10^5;
end

nbScl = length(data);

for s = 1:nbScl
    % data_noisy = P.add_white_noise_global(data{m,ss}, nlvl);
    data_noisy = obj.add_white_noise_global(data{s}, nlvl);

    out = obj.reconstruct_CGPT(data_noisy.MSR_noisy, ord, maxiter, tol, symmode, method, op);
    
    CGPTt{s} = tools.cell2mat3D(out.CGPT);
end
