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
    As = PDE.Conductivity_R2.make_matrix_WPT(src, cfg.center, WF, nbScl, stdmode);
    [Asn, Ws] = tools.col_normalization(As);
    
    % The matrix related to receivers
    if cfg.fixed_rcv
        % If the receiver does not depend on the source, use a simplified model
        % for speed
        if cfg.coincided_src_rcv
            Ar = As;
        else
            Ar = PDE.Conductivity_R2.make_matrix_WPT(cfg.rcv(1), cfg.center, WF, nbScl, stdmode);
        end
        [Arn, Wr] = tools.col_normalization(Ar);
        
        L = @(x,tflag)tools.linsys.SXR_op(x,As,Ar,tflag); % linear operator
        Ln = @(x,tflag)tools.linsys.SXR_op(x,Asn,Arn,tflag); % normalized linear operator
    else
        for n=1:cfg.Ns
            Ar{n} = PDE.Conductivity_R2.make_matrix_WPT(cfg.rcv(n), cfg.center, WF, nbScl, stdmode);
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
