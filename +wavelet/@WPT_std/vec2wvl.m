%% TODO
function WPT = vec2wvl(nbDir, nbScl,  X0)
% Convert a WPT vector to a WPT object.

if ~ismatrix(X0)
    error('The input argument must be a matrix or vector!');
end

N = sqrt(numel(X0));
X = reshape(X0, N, N);

% XAD = X(end-nbApprx+1:end,1:end-nbApprx);
% XDA = X(1:end-nbApprx,end-nbApprx+1:end);
% XDD = X(1:end-nbApprx,1:end-nbApprx);

DD = cell(nbDir);
DA = cell(nbDir,1);
AD = cell(1,nbDir);
% WPT.Coeff = cell(nbDir*nbScl+1);

for lr=1:nbDir
    for lc=1:nbDir
        toto = X((lr-1)*N+1:lr*N, (lc-1)*N+1:lc*N);
        
        % Loops on the scale
        rr=0;
        %        for j1=Jmin:Jmax
        % jr = j1-Jmin + 1;
        for jr=1:nbScl
            nr = prod(dmesh.shape{jr});
            cc=0;
            
            for jc=1:nbScl
                %                jc = j2-Jmin+1;
                nc = prod(dmesh.shape{jc});
                
                blk = toto(rr+1:rr+nr, cc+1:cc+nc);
                DD{lr,lc}{jr,jc} = blk;
                %                Coeff{(lr-1)*nbScl+jr,(lc-1)*nbScl+jc} = blk;
                
                cc = cc+nc;
            end
            rr = rr+nr;
        end
    end
end

% Detail-Apprx
%
% Loops on the direction
for lr=1:nbDir
    toto = X((lr-1)*N+1:lr*N, end-nbApprx+1:end);
    rr=0;
    
    % Loops on the scale
    for j=Jmin:Jmax
        jr = j-Jmin + 1;
        nr = prod(dmesh.shape{jr});
        
        blk = toto(rr+1:rr+nr, :);
        DA{lr}{jr} = blk;
        Coeff{(lr-1)*nbScl+jr,end} = blk;
        
        rr = rr+nr;
    end
end

% Apprx-Detail
%
% Loops on the direction
for lc=1:nbDir
    toto = X(end-nbApprx+1:end,(lc-1)*N+1:lc*N);
    cc=0;
    
    % Loops on the scale
    for j=Jmin:Jmax
        jc = j-Jmin + 1;
        nc = prod(dmesh.shape{jc});
        
        blk = toto(:, cc+1:cc+nc);
        AD{lc}{jc} = blk;
        Coeff{end, (lc-1)*nbScl+jc} = blk;
        
        cc = cc+nc;
    end
end

AA = X(end-nbApprx+1:end,end-nbApprx+1:end);
Coeff{end,end} = AA;
end