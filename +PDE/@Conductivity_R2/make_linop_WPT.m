function out = make_linop_WPT(cfg, WF, nbScl, stdmode)
% out = make_linop_WPT(cfg, WF, nbScl, stdmode);
% Make linear operator of the WPT system
% INPUTS:
% cfg: configuration of acquisition system, an object of acq.mconfig class
% WF: object of wavelet.OrthoWvl
% nbScl: number of if nbScl=0 then only approximation coefficients are used
% stdmode: true for standard repersentation
% OUTPUT: a structure including the following fields
% As, Ar: source/receiver matrices
% L: linear operator L(X) = As*X*Ar'
% Ws,Wr: column weight of As and Ar
% W: the diagonal weight such that L(X) = Ln(W.*X)
% Asn, Arn: column normalized version of As and Ar
% Ln: column normalized version of L
%
% Remark:
% For the linear operator L, its matrix form has the (p,q)-th column vector
% (As_{mp}*Ar_{nq})_{mn}, therefore, the column-normalized operator is L(X) = Bs*X*Br', with Bs,
% Br the column normalized matrices of As, Ar.
    
    if nargin < 4
        stdmode = 1;
    end
    if nargin < 3
        nbScl = 0;
    end

    % The matrix related to sources:
    src = cfg.all_src();
    As = make_matrix_WPT(src, cfg.center, WF, nbScl, stdmode);
    [Asn, Ws] = tools.col_normalization(As);
    
    % The matrix related to receivers
    if cfg.fixed_rcv
        % If the receiver does not depend on the source, use a simplified model
        % for speed
        if cfg.coincided_src_rcv
            Ar = As;
        else
            Ar = make_matrix_WPT(cfg.rcv(1), cfg.center, WF, nbScl, stdmode);
        end
        [Arn, Wr] = tools.col_normalization(Ar);
        
        L = @(x,tflag)tools.linsys.SXR_op(x,As,Ar,tflag); % linear operator
        Ln = @(x,tflag)tools.linsys.SXR_op(x,Asn,Arn,tflag); % normalized linear operator
    else
        for n=1:cfg.Ns
            Ar{n} = make_matrix_WPT(cfg.rcv(n), cfg.center, WF, nbScl, stdmode);
            [Arn{n}, Wr{n}] = tools.col_normalization(Ar{n});            
        end
        L = @(x,tflag)tools.linsys.SXR_op_list(x,As,Ar,tflag); % linear operator
        Ln = @(x,tflag)tools.linsys.SXR_op_list(x,Asn,Arn,tflag); % normalized linear operator
    end

    out.L = L;
    out.As = As;
    out.Ar = Ar;

    out.Ln = Ln;
    out.Asn = Asn;
    out.Arn = Arn;
    
    out.Ws = Ws;
    out.Wr = Wr;
    
    if iscell(Wr)
        out.W = sqrt(abs(Ws').^2 * abs(cell2mat(Wr)).^2);
    else
        out.W = reshape(Ws, [], 1) * reshape(Wr, 1, []);
    end

    % Lop.times = @(x) L(x,'notransp');
    % Lop.trans = @(x) L(x,'transp');
    % out.Lop = Lop;

    % Lnop.times = @(x) Ln(x,'notransp');
    % Lnop.trans = @(x) Ln(x,'transp');
    % out.Lnop = Lnop;
end

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
