classdef WPT_std
% Class for wavelet Polarization Tensor (WPT) of standard representatieon

    properties(SetAccess=protected)
        nbScl % number of scale of detail spaces
        nbDir % number of different directional wavelets in each scale of detail spaces

        DD % Detail-Detail coefficients, 2D cell of dimension nbScl-by-nbScl, organized as:
           %
           % DD{k1,k2}{j1,j2} for k1,k2=1..nbDir and j1,j2=1..nbScl, is a N1-by-N2
           % matrix, with the {m,n}-th entry the WPT between the wavelet
           % psi^{k1}_{j1,m} and the wavelet psi^{k2}_{j2,n}, and N1 is the number of
           % active wavelets in the scale {j1,k1} (similar for N2).

        DA % Detail-Approximation coefficients, 1D cell of dimension nbScl-by-1, organized
           % as:
           %
           % DA{k1}{j1} for k1=1..nbDir and j1=1..nbScl is a N1-by-N2 matrix, with
           % the {m,n}-th entry the WPT
           % between the wavelet psi^{k1}_{j1,m} and the scaling function phi_{Jmax,n},
           % and N1 is the number of active wavelets in the scale {j1,k1}, and N2
           % is the the number of active scaling functions in the scale Jmax.

        AD % Approximation-Detail coefficients, 1D cell of dimension 1-by-nbScl, organized
           % as:
           %
           % AD{k1}{j1} for k1=1..nbDir and j1=1..nbScl is a N1-by-N2 matrix, with
           % the {m,n}-th entry the WPT between the scaling function phi_{Jmax,m} and
           % the wavelet psi^{k1}_{j1,n}, and N1 is the the number of active scaling
           % functions in the scale Jmax, and N2 is the number of active wavelets in
           % the scale {j1,k1}.

        AA % Approximation-Approximation coefficients, N1-by-N1 matrix, with
           % the {m,n}-th entry the WPT between the scaling function phi_{Jmax,m}
           % and phi_{Jmax,n}, and N1 is the the number of active scaling
           % functions in the scale Jmax.
        
        datac % whole coefficient matrix in cell format        
              % Format of datac:
              % Block row 1:
              %
              % [ W^{1,1}_{Jmin,Jmin}...W^{1,nbDir}_{Jmin,Jmin},
              % W^{1,1}_{Jmin,Jmin+1}...W^{1,nbDir}_{Jmin,Jmin+1}...
              % W^{1,1}_{Jmin,Jmax}...W^{1,nbDir}_{Jmin,Jmax}, W^{1,0}_{Jmin,Jmax} ] 
              %
              % ...  
              % Block row nbDir:
              %
              % [ W^{nbDir,1}_{Jmin,Jmin}...W^{nbDir,nbDir}_{Jmin,Jmin},
              % W^{nbDir,1}_{Jmin,Jmin+1}...W^{nbDir,nbDir}_{Jmin,Jmin+1}...
              % W^{nbDir,1}_{Jmin,Jmax}...W^{nbDir,nbDir}_{Jmin,Jmax},
              % W^{nbDir,0}_{Jmin,Jmax} ]
              %
              % Block row nbDir+1:
              %
              % [ W^{1,1}_{Jmin+1,Jmin}...W^{1,nbDir}_{Jmin+1,Jmin},
              % W^{1,1}_{Jmin+1,Jmin+1}...W^{1,nbDir}_{Jmin+1,Jmin+1}...
              % W^{1,1}_{Jmin+1,Jmax}...W^{1,nbDir}_{Jmin+1,Jmax}, W^{1,0}_{Jmin+1,Jmax} ]
              %
              % etc.
              % 
              % W^{k1,k2} is DD{k1,k2} if k1,k2=1..nbDir, DA{k1} if k2=0, AD{k2} if k1=0, and AA if k1=k2=0.
        
        data % matrix version of datac
    end
    
    methods
        function obj = WPT_std(DD,DA,AD,AA)
            obj.DD = DD;
            obj.DA = DA;
            obj.AD = AD;
            obj.AA = AA;
            obj.nbDir = size(DD, 1);
            if isempty(DD)
                obj.nbScl = 0;
            else
                obj.nbScl = size(DD{1,1}, 1);
            end
        end
        
        function val = lowscls(obj, nbScl)
        % Keep only coarse scale coefficients
            if obj.nbDir>0 && obj.nbScl>0
                DD = cell(obj.nbDir); 
                DA = cell(obj.nbDir,1); 
                AD = cell(1,obj.nbDir); 
                
                j0 = obj.nbScl-nbScl+1;
                
                for kr=1:obj.nbDir
                    DA{kr} = obj.DA{kr}(j0:end, :);
                    AD{kr} = obj.AD{kr}(:, j0:end);
                    for kc=1:obj.nbDir                
                        DD{kr,kc} = obj.DD{kr,kc}(j0:end, j0:end);
                    end
                end
                
                val = wavelet.WPT_std(DD,DA,AD,obj.AA);
            else
                val = wavelet.WPT_std({},{},{},obj.AA);
            end
        end

        function val = onedir(obj, kr, kc)
        % Keep only coefficients of one direction
            if obj.nbScl == 0
                val = wavelet.WPT_std([], [], [], obj.AA);
            else
                DD = obj.DD(kr,kc);
                DA = obj.DA(kr);
                AD = obj.AD(kc);
                val = wavelet.WPT_std(DD,DA,AD,obj.AA);
            end        
        end
        
        function M = mDD(obj, kr, kc)
            M = cell2mat(obj.DD{kr,kc});
        end
        
        function M = mDA(obj, kr)
            M = cell2mat(obj.DA{kr});
        end
        
        function M = mAD(obj, kc)
            M = cell2mat(obj.AD{kc});
        end
        
        function C = get.data(obj)
            C = cell2mat(obj.datac);
        end
        
        function C = get.datac(obj)
            C = cell(obj.nbDir*obj.nbScl+1);
            for lr=1:obj.nbDir
                for lc=1:obj.nbDir
                    C(lr:obj.nbDir:end-1, lc:obj.nbDir:end-1) = obj.DD{lr,lc};                
                end 
            end
            
            for lr=1:obj.nbDir
                C(lr:obj.nbDir:end-1, end) = obj.DA{lr};
            end
            
            for lc=1:obj.nbDir
                C(end, lc:obj.nbDir:end-1) = obj.AD{lc};
            end

            C{end,end} = obj.AA;
            
            % for lr=1:obj.nbDir
            %     for lc=1:obj.nbDir
            %         % Loops on the scale
            %         for Jr=Jmin:Jmax
            %             jr = Jr-Jmin + 1;   
            %             %
            %             for Jc=Jmin:Jmax
            %                 jc = Jc-Jmin+1;
                            
            %                 C{(jr-1)*obj.nbDir+lr, (jc-1)*obj.nbDir+lc} = obj.DD{lr, ...
            %                                     lc}{jr,jc};
            %             end
            %         end
            %     end
            % end
            
            % for lr=1:obj.nbDir
            %     % Loops on the scale    
            %     for Jr=Jmin:Jmax
            %         jr = Jr-Jmin + 1;
            %         C{(jr-1)*obj.nbDir+lr,end} = obj.DA{lr}{jr};
            %     end
            % end
            
            % for lc=1:obj.nbDir
            %     % Loops on the scale
            %     for Jc=Jmin:Jmax
            %         jc = Jc-Jmin + 1;
            %         C{end, (jc-1)*obj.nbDir+lc} = obj.AD{lc}{jc};
            %     end
            % end

            % C{end,end} = obj.AA;
        end
    
    end
    
    methods(Static)
        WPT = vec2wvl(X)
        
        function X = diag_imaging(W)
        % X = diag_imaging(W)
        % Imaging by taking the diagonal of WPT. 
        % INPUT:
        % W: WPT coefficient matrix, of dimension N^2 X N^2
        % OUTPUT:
        % X: an image of dimension N X N.
            
            N = sqrt(numel(W));
            W = reshape(W, N, N);
            M = round(sqrt(N));

            if M^2 ~= N
                error('Dimension error!');
            end
            
            X = abs(full(reshape(diag(W), M, M)));
        end
        
        function X = max_imaging(W, nw)
        % X = max_imaging(W, nw)
        % Imaging by maxima of WPT. 
        % INPUT:
        % W: WPT coefficient matrix, of dimension N^2 X N^2
        % nw: use the first n maxima for imaging
        % OUTPUT:
        % X: an image of dimension N X N.

            if nargin<2
                nw=1;
            end
            
            N = sqrt(numel(W));
            W = reshape(W, N, N);
            M = round(sqrt(N));

            if M^2 ~= N
                error('Dimension error!');
            end
            
            X = zeros(M^2,1);

            for m=1:M^2
                % Recall the definition of W matrix:
                % W(m,n) = \int_{\p D} psi_n (\lambda I-K_D^*)^-1 [d psi_m dn] ds(x)

                if nw == 1
                    [v, idx] = max(abs(W(m,:))); % W(m,:) takes the m-th source

                    X(idx) = X(idx) + (abs(W(m,idx)));
                    % X(idx) = X(idx)+v;
                    % X(idx) = max(X(idx), v);                    
                else
                    [v, idx] = sort(abs(W(m,:)), 'descend');

                    X(idx(1:nw)) = X(idx(1:nw)) + mean(v(1:nw));
                    % X(idx(1:nw)) = X(idx(1:nw)) + mean(abs(W(:,idx(1:nw))));
                    % X(idx(1:nw)) = X(idx(1:nw))+v(1:nw)';                    
                    % X(idx(1:nw)) = max(X(idx(1:nw)), v(1:nw)');                    
                end
            end
            X = full(reshape(X, M, M));
        end
        
        M = make_WPT_mask(N, n)
    end        
end

