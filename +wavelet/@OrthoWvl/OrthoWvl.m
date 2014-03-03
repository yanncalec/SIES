classdef OrthoWvl < wavelet.Frame
% 2D Orthogonal wavelet class with periodic boundary condition
% 
% We introduce an augmented ROI (AROI) which is one scale finer than the scale Jmin (Jnum=Jmin-1),
% and slightly larger than the ROI so that wavelets of scales from Jmin to Jmax (as well as scaling
% functions of scale Jmax) intersecting with ROI do not intersect the boundary of AROI.
%
% For a function f defined on the AROI, its point values (up to a constant factor) at the grid of
% the scale Jnum are used as input to the wavelet decomposition, and the fast transform algorithm
% decomposes it til the scale Jmax. At each scale the decomposition increases the dimension (due to
% the non-periodic convolution, the transform it is not orthogonal). At the synthesis the dimension
% will also increase, but the reconstruction limited on the grid of AROI is identical to the
% original f. What is synthesized with the coeffcients of scales from Jmax to Jmin, is the
% approximation coeffcients of the scale Jmin-1=Jnum.
% 
% Furthermore, the position number n\in\Z^2 of wavelet coefficients are tracked during the
% decomposition. At each scale the ranges of the position index are stored in wrangex and wrangey,
% which are necessary when we do interpolation with wavelets. The wavelet contributing to the ROI
% (called active) are coded as a 0-1 mask and stored in wmask. wmask.A has the same dimension as the
% approximation coefficient (A returned by analysis), similarly wmask.D{j,k} has the same dimension
% as detail coefficients (scale index j, direction k=1,2,3). To make it clear: wmask has the same
% dimension as that of the wavelet coefficient, and ApprxSpace{n} records the dimension of the
% active wavelet of a scale (a submatrix of wmask.A or wmask.D)
%
% Convention: A 2D image is the value of a function F evaluated on a 2d grid (e.g. inside
% [-1,1]X[-1,1]). The grid is generated using the function 'meshgrid' by increasing in both x and y
% axis, e.g. [Sx, Sy]=meshgrid(-1:0.1:1, -1:0.1:1). Visually, the 2D image is flipped up-down. The
% same order of representation must be respected in all calculations, since the wavelet transform of
% a flipped up-down or left-right image is not the same (not even in reversed order) of the original
% one. It's only at the last visulization step that we correct the image using the function
% flipud. See for example ShiftInvSpace.mesh for details.
%

    properties(SetAccess=protected)
        % quadrature mirror filter
        hcoeff % low pass filter
        hrange % index range of the low pass filter
        gcoeff % high pass filter
        grange % index range of the high pass filter

        Jnum % scale number of numerical computation, Jnum<=Jmin-1. Jnum defines
             % the finest grid for approximation of a L2(R^2) function.
        
        nbDir = 3 % number of different types of wavelets (directions) in each detail space

        wvlord % order of the wavelet, eg, for Daubechie wavelet this can be 1,2,...

        Tphi % table of values of 1D mother scaling function
        Tpsi % table of values of 1D mother wavelet function
        dTphi % table of values of the derivative of 1D mother scaling function
        dTpsi % table of values of the derivative of 1D mother wavelet function
        Gphi % regular grid on which phi is evaluated
        Gpsi % regular grid on which psi is evaluated
        
        phisupp % interval of support of the 1D scaling function phi
        psisupp % interval of support of the 1D wavelet function psi

        wvlsize % support size of the 1D mother wavelet psi or scaling function phi
        
        wvlsize_Jmin % support size of the wavelet at the finest scale
        wvlsize_Jmax % support size of the wavelet at the coarsest scale

        AROI % approximation space for augmented ROI corresponding to the scale Jnum
             %
             % Instance of wavelet.ROI class
        
        wmask % mask for active wavelet coeffcients
              %
              % wmask.A: mask for the approximation scale
              % wmask.D{j,k} mask for the detail space at the scale j and direction k
              % wmask.all: concatenation of [A, D{1,1}, D{1,2}, D{1,3}...] in a vector

        wrangex % position range in x-axis of wavelet coefficients
        wrangey % position range in y-axis of wavelet coefficients
        
        wptr % starting index pointer for detail wavelet coefficients
             %
             % The wavelet coefficient vector is the concatenation W=[A, D{1,1},
             %D{1,2}, D{1,3}, D{2,1}...]. So wptr(1) is the index where D{1,1} begins
             %in W, and wptr(4) is for D{2,1}, etc.        
    end
    
    properties(Access=protected)
        ordmax % maximum monomial order used in numerical approximation of inner product <f,psi>
        
        monomcoeff % all coefficients monomcoeffs{k}(a1,a2) = <x^a, psi^k> 
                   % for a=(a1,a2), a1,a2>=wvlord and a1,a2<=ordmax,
                   % and k=1,2,3
        
        wgIter = 12 % number of iteration in the algorithm for computing wavelet graph

        JnumPow = 11 % in the finest approximation space, the ROI is sampled by
                     % 2^JnumPow X 2^JnumPow points
    end
    
    methods        
        function obj = OrthoWvl(wvlfamily, wvlord, ROI_center, ROI_width, Rsl)
        % INPUTS:
        % wvlfamily: type of wavelet family, 'Daubechies' for Daubechies wavelet. See wavelet.MakeONFilter
        % wvlord: order of the wavelet. See wavelet.MakeONFilter
        % ROI_center: center of the square ROI
        % ROI_width: width of the ROI
        % Rsl: resolution of the imaging system
            
            obj.Rsl = Rsl; 
            obj.sclFct = 2; % Scaling factor is 2 for orthogonal wavelets
            obj.splStep = 1; % No over-sampling

            obj.wvlname = [wvlfamily, '_', num2str(wvlord)];
            obj.wvlord = wvlord; 
            
            obj.hcoeff = wavelet.MakeONFilter(wvlfamily, wvlord); % Get the qmf
            obj.hrange = [0, length(obj.hcoeff)-1]; % Range of the qmf. This range is valid for Daubechies wavelet.
            
            % Calculate the wavelet's graph. The number of iteration should not be smaller than
            % 11, otherwise the direct computation of <f, Psi> by quadrature would be inaccurate.
            [obj.Tphi, obj.Tpsi, obj.Gphi, obj.Gpsi, obj.dTphi, obj.dTpsi] = ...
                wavelet.OrthoWvl.wvlgraph(obj.hcoeff, obj.gcoeff, obj.wgIter); 

            % Determine the support size of the given wavelet
            obj.wvlsize = obj.Gphi(end)-obj.Gphi(1);
            
            % The coarsest scale is determined from the system's resolution and
            % the ROI size:
            obj.Jmax = min(fix(log2(obj.Rsl)), fix(log2(ROI_width/obj.wvlsize)));
            % obj.Jmax = max(fix(log2(obj.Rsl)), fix(log2(ROI_width/obj.wvlsize)));

            % Augmented ROI (AROI) is determined so that any wavelet of the scale Jmax
            % intersecting the border of AROI does not meet the ROI. In this way
            % we avoid the bord effet on ROI.
            awidth = 2*(ROI_width/2 + 1.1*obj.wvlsize_Jmax); % width

            % Scale for numerical computation
            obj.Jnum = min(obj.Jmax-1, fix(log2(ROI_width/2^obj.JnumPow))); %
            % obj.Jnum = min(obj.Jmax-1, -8);

            % Shift invariant space for approximation on AROI of scale Jnum
            obj.AROI = wavelet.ShiftInvSpace(ROI_center, awidth, awidth, obj.Jnum, 0, ...
                                             obj.hrange);
            
            % The finest scale is one scale above Jnum (coarser than Jnum)
            obj.Jmin = obj.Jnum+1;

            % Make shift-invarince spaces
            obj.DetailSpace = cell(obj.nbScl, 3);
            for j = 1:obj.nbScl
                J = obj.Jmax-j+1;
                obj.ApprxSpace{j} = wavelet.ShiftInvSpace(ROI_center, ROI_width, ROI_width, J, 0, ...
                                                          obj.hrange);
                for k=1:3                    
                    obj.DetailSpace{j,k} = wavelet.ShiftInvSpace(ROI_center, ROI_width, ROI_width, ...
                                                                 J, k, obj.hrange);
                end
            end
            % Add one more approximation space, corresponding to the scale Jnum (in the same
            % scale as AROI, but will be used only for defining active wavelets).
            obj.ApprxSpace{obj.nbScl+1} = wavelet.ShiftInvSpace(ROI_center, ROI_width, ROI_width, ...
                                                              obj.Jnum, 0, obj.hrange);
            % the ROI is actually also the finest approximation space
            obj.ROI = obj.ApprxSpace{obj.nbScl+1}; 

            % Make masks for active wavelet coeffcients
            X0 = ones(obj.AROI.dim);
            [A, D] = obj.analysis(X0); % Perform a test decomposition
            
            % The position range of approximation coeffs is larger than that of active
            % wavelets. By extending the range we make a mask for active coefficients:
            obj.wmask.A = obj.ApprxSpace{1}.extend_mask(A.rangex, A.rangey);
            
            % Range and dimension of the approximation coefficients
            obj.wrangex.A = A.rangex; obj.wrangey.A = A.rangey;
            obj.wmask.all = [obj.wmask.A(:)];
            wptr = [1 1+numel(obj.wmask.A)];

            % Do the same thing for detail spaces
            for j=1:obj.nbScl
                for k=1:obj.nbDir
                    obj.wmask.D{j,k} = obj.DetailSpace{j,k}.extend_mask(D{j,k}.rangex, ...
                                                                      D{j,k}.rangey);
                    obj.wrangex.D{j,k} = D{j,k}.rangex;
                    obj.wrangey.D{j,k} = D{j,k}.rangey;
                    wptr = [wptr wptr(end)+numel(obj.wmask.D{j,k})];
                end
                obj.wptr = wptr(2:end-1);
                obj.wmask.all = [obj.wmask.all; obj.wmask.D{j,1}(:); obj.wmask.D{j,2}(:); obj.wmask.D{j,3}(:)];
            end

            % % Compute wavelet coefficients for monomials:
            % obj.ordmax = obj.wvlord+5;
            % obj.monomcoeff{1} = wavelet.OrthoWvl.monomcoeff_mex(obj.Tphi, obj.phisupp, obj.Tpsi, ...
            %                                                   obj.psisupp, obj.wvlord, obj.ordmax);
            % obj.monomcoeff{2} = wavelet.OrthoWvl.monomcoeff_mex(obj.Tpsi, obj.psisupp, obj.Tphi, ...
            %                                                   obj.phisupp, obj.wvlord, obj.ordmax);
            % obj.monomcoeff{3} = wavelet.OrthoWvl.monomcoeff_mex(obj.Tpsi, obj.psisupp, obj.Tpsi, ...
            %                                                   obj.psisupp, obj.wvlord, obj.ordmax);
            
        end        

        function val = get.grange(obj)
            val = 1 - [obj.hrange(2), obj.hrange(1)]; % the mathematical definition
        end
        
        function val = get.gcoeff(obj)
            val = wavelet.OrthoWvl.qmf_h2g(obj.hcoeff, obj.hrange);
        end
        
        function val = get.phisupp(obj)
            val=[obj.Gphi(1), obj.Gphi(end)];
        end
        
        function val = get.psisupp(obj)
            val=[obj.Gpsi(1), obj.Gpsi(end)];
        end        
        
        function val = get.wvlsize_Jmin(obj)
            val = obj.wvlsize * 2^obj.Jmin;
        end
        
        function val = get.wvlsize_Jmax(obj)
            val = obj.wvlsize * 2^obj.Jmax;
        end

        function [rx, ry] = getsupport(obj, J, k, n)
        % Get the support of a wavelet psi^k_{J,n}
            n=reshape(n, 1, 2);
            if k==0
                rx = 2^J*(obj.phisupp + n);
                ry = 2^J*(obj.phisupp + n);
            elseif k==1
                rx = 2^J*(obj.phisupp + n);
                ry = 2^J*(obj.psisupp + n);
            elseif k==2
                rx = 2^J*(obj.psisupp + n);
                ry = 2^J*(obj.phisupp + n);
            else
                rx = 2^J*(obj.psisupp + n);
                ry = 2^J*(obj.psisupp + n);
            end            
        end
    end
    
    %% Abstract methods
    methods
        function [A, D, W] = analysis(obj, X)
        % [A, D, W] = analysis(obj, X)
        %
        % Wavelet decomposition applied on AROI. The transform is performed from the scale Jnum
        % (AROI grid) to Jmax.
        %
        % INPUTS:
        % X: value of a function on the AROI grid (or the approximation coefficient of scale Jnum), 2D
        % matrix of dimension AROI.dim
        %
        % OUTPUTS:
        % A,D: two structures
        % A.coeff: approximation coefficient
        % A.rangex, A.rangey: position range in x and y-axis of A.coeff
        % D{j,k}.coeff, D{j,k}.rangex(y): same as A but for detail coeffs
            
            [nrow, ncol] = size(X);
            if nrow~=obj.AROI.dim(1) || ncol~=obj.AROI.dim(2)
                error('Dimension error: input matrix must be have the same dimension as the AROI grid');
            end
            
            Nx1 = obj.AROI.rangex(1);
            Ny1 = obj.AROI.rangey(1);

            % number of the decompisition level, remark that nlevel == obj.nbScl
            nlevel = obj.Jmax - obj.Jnum; 
            [A, D] = wavelet.OrthoWvl.wavedec2(X, Nx1, Ny1, obj.hcoeff, obj.hrange, ...
                                               nlevel);
            W = wavelet.OrthoWvl.scl2vec(A, D);
        end
        
        function X = synthesis(obj, A, D)
        % X = synthesis(obj, A, D) 
        % Synthesize a function on the AROI from its wavelet coefficients A,D returned by analysis(X). The
        % output X has the same dimension as the AROI grid.
        % 
            A0 = wavelet.OrthoWvl.waverec2(A, D, obj.hcoeff, obj.hrange);
            % The size of reconstructed signal can be much larger than the
            % original one due to the convolution (with many extra zeros on
            % borders for example)
            X = obj.AROI.resize(A0.coeff, A0.rangex, A0.rangey, 'crop');
        end
        
        function [A, D, W] = analysis_ROI(obj, X)
        % [A, D, W] = analysis_ROI(obj, X)
        %
        % Wavelet decomposition restricted on the ROI. Inactive coefficients (those not contributing to ROI)
        % are set to zero.
            
            [A, D, W] = obj.analysis(X);
            W = W .* obj.wmask.all; % keep active wavelets
            [A, D] = obj.vec2scl(W);
        end
        
        function X = synthesis_ROI(obj, A, D)
        % X = synthesis_ROI(obj, A, D)
        %
        % Synthesize a function on a ROI from its wavelet coefficients A,D. The scale of the shift invariant
        % space of approximation is automatically determined from the number of decomposition levels in D,
        % and the output X has the same dimension as the ROI grid.
        % 
            if size(D,1)>obj.nbScl
                error('Wrong number of levels in D.');
            end
            
            A0 = wavelet.OrthoWvl.waverec2(A, D, obj.hcoeff, obj.hrange);

            % size(D,1) is the number of decomposition level, obj.Jmax-size(D,1) is the scale of synthesis.
            j0 = size(D,1)+1; % index of the scale of synthesis
            X = obj.ApprxSpace{j0}.resize(A0.coeff, A0.rangex, A0.rangey, 'crop');
        end
        
        function [A,D] = dual_analysis(obj, X)
            [A,D] = obj.analysis(X);
        end

        function X = dual_synthesis(obj, A, D)
            X = obj.synthesis(A,D);
        end            
        
        val = eval_apprx(obj, X)
        
        val = eval_detail(obj, X, J, l)
        
        % Resynthesis of a function on the ROI by using coefficients of coarse scales.
        [Xr, err] = apprx_ROI(obj, W, nbScl, Xim)

        %% Utility functions
        
        [A, D] = vec2scl(obj, W) % Convert a wavelet coefficient vector to [A,D] form        
        W = invert_scl(obj, A, D) % Invert the scale order in wavelet coefficients

        function Winv = invert_scl_v(obj, W, nbScl)
        % Interface function to invert_scl
            [A, D] = obj.vec2scl(W);
            Winv = obj.invert_scl(A, D(1:nbScl,:));
        end
        
        % Mask for extraction of wavelet coefficients (active of not) of a given scale and direction
        mask = mask_onescl(obj, j, k)        
        % Mask for extraction of wavelet coefficients (active of not) up to a given scale
        mask = mask_lowscls(obj, nbScl)                
        % Mask for extraction of active wavelet coefficients of a given scale and direction
        mask = mask_active_onescl(obj, j, k)
        % Mask for extraction of active wavelet coefficients up to a given scale        
        mask = mask_active_lowscls(obj, nbScl)            
        
        % Extact wavelet coefficients (active or not) of the scale (j,k)
        val = extract_coeffs(obj, W, j, k) % output is a matrix
        val = extract_coeffs_full(obj, W, j, k) % this preserves the dimension of W
        
        function Wa = extract_coeffs_lowscls(obj, W, nbScl)
        % Extract the first nbScl scales from a coefficient vector

            mask = obj.mask_lowscls(nbScl);
            Wa = W .* mask;
        end            

        % Extact only activ wavelet coefficients of the scale (j,k)
        function val = extract_active_coeffs(obj, W, j, k)
            if nargin < 4
                if j ~= 0
                    error(['The direction index k must be given for detail ' ...
                           'spaces']);
                end
                k = 0;
            end

            mask = obj.mask_active_onescl(j,k);
            val = W(find(mask));
        end
            
        function val = extract_active_coeffs_lowscls(obj, W, nbScl)
            mask = obj.mask_active_lowscls(nbScl);
            val = W(find(mask));
        end

        % Make wavelet interpolation matrix on a ROI
        [Phi, dPhiX, dPhiY] = interpolation_mat_ROI(obj, X, J) 
        
        function W = WPT_std(obj, D, kappa, nbScl)            
        % admittivity: kappa = cnd + 1i * pmtt * freq
            
            W = wavelet.OrthoWvl.theoretical_WPT_std(D.points, D.tvec, D.normal, D.avec, D.sigma, ...
                                                     kappa, D.center_of_mass, obj.Jmax, ...
                                                     obj.ApprxSpace{1}, obj.DetailSpace(1:nbScl,:), ...
                                                     obj.Tpsi, obj.Tphi, obj.Gpsi, obj.Gphi);
        end
        
        function W = WPT_nonstd(obj, D, freq)
            W=[];
        end        
        
        %% For Green function of Laplacian
        [Wa, Xim] = decomp_Green2D(obj, Z) % Decomposition the Green function of Laplacian by
                                           % fast wavelet transform

        [Wa, Xim] = decomp_Green2D_num(obj, Z, nbScl) % Same as decomp_Green2D but by direct numerical integration

        A = Green2D_apprxcoeff(obj, Z)
        D = Green2D_detailcoeff(obj, Z, nbScl)
        
        % function Wa = Green2D_apprxcoeff(obj, Z)
        %     X = tools.Laplacian.Green2D(obj.AROI.meshpoints, Z(:)) * 2^obj.Jnum;
        %     % decomposition from Jnum to Jmax
        %     [~, ~, W] = obj.analysis(reshape(X, obj.AROI.dim)); 
        %     % W = W(1:numel(obj.wdim.A));
        %     Wa = W(find(obj.wmask.A(:)));
        % end
        
        % Put the coefficient vector into the classical wavelet tree matrix form
        M = coeff2mat(obj, C) 
        
        W = WPT_vec2wvl(obj, X) % convert a vector to a WPT structure
        X = WPT_wvl2vec(obj, W) % convert a WPT structure to a vector
    end

    methods(Static)        
        function [rangex, rangey] = matrix_idx_range(C, Nx1, Ny1)
        % [rangex, rangey] = matrix_idx_range(C, Nx1, Ny1)
        % Compute the physical range index for a matrix C, whose the physical row
        % index starting by Ny1, and physical column starting by Nx1.
            
            [nrow, ncol] = size(C);
            rangex = [Nx1, Nx1+ncol-1]; % column
            rangey = [Ny1, Ny1+nrow-1]; % row
            
            % [rangex, rangey] = meshgrid(Nx1:Nx1+ncol-1, Ny1:Ny1+nrow-1);
        end

        function [gcoeff, grange] = qmf_h2g(hcoeff, hrange)
        % [gcoeff, grange] = qmf_h2g(hcoeff, hrange)
        % Compute the conjugated mirror filter g from a filter h.
        % hcoeff: the filter h
        % hrange: the starting and ending index of h
            grange = 1 - [hrange(2), hrange(1)];
            sg = (-1).^(hrange(1):hrange(2));
            gcoeff = fliplr(conj(hcoeff) .* sg); 

            % sg = (-1).^(1-(grange(1):grange(2)));
            % gcoeff = conj(hcoeff(end:-1:1)) .* sg;             
        end
        
        W = scl2vec(A, D, nbScl) % convert wavelet coefficients [A,D] to a vector form
        
        [V, dV] = lookup_table_wvl_mex(X, N, Gx1, Gx2, dx, Tx, dTx)

        A = Green2D_apprxcoeff_mex(Tphi, phisupp, Jmax, rangex, rangey, Z)
        D = Green2D_detailcoeff_mex(Tphi, phisupp, Tpsi, psisupp, J, rangex, rangey, Z)
        C = monomcoeff_mex(Tpsi, psisupp, Tpsi, psisupp, wvlord, ordmax)
                
        [Tphi, Tpsi, Gphi, Gpsi, dTphi, dTpsi] = wvlgraph(hcoeff, gcoeff, niter)
        [Nx, Ny] = active_position(Z0, Lx, Ly, N1, N2, J, k)
        
        [A, D] = wavedec2(X, Nx1, Ny1, hcoeff, hrange, nbScl)
        A0 = waverec2(A, D, hcoeff, hrange)

        [A, D1, D2, D3] = dwt2(X, Nx1, Ny1, hcoeff, hrange, gcoeff)
        out = idwt2(A, D1, D2, D3, hcoeff, hrange, gcoeff)
        
        [val, dX, dY] = evaluation(X, J, N, Tx, Ty, Gx, Gy)
        
        WPT = theoretical_WPT_std(D, tvec, normal, avec, sigma, kappa, Z0, Jmax, ...
                                  ApprxSpace, DetailSpace, Tpsi, Tphi, Gpsi, Gphi)
    end

end
