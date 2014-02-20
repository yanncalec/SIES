classdef ShiftInvSpace < wavelet.ROI
% Class for shift invariant space 
%
% Mathematically speaking, a shift invariant space is the linear space generated by translating a
% mother function, and an ideal design would include the evaluation of the mother
% function. Notheless, we decide to put the evaluation part into the sub class of |wavelet.Frame|,
% and here only includes properties like positions of the mother function (which is of course a
% awkward design).
%   
    
    properties(SetAccess=protected)
        J % scale number
        id % identification: 0 for an approximation space, 1..k for directional detail spaces

        ortho % true for orthogonal wavelet basis, false for redundant basis
        qmf_range = [] % index range [N1, N2] of the quadratic mirror filter
        
        rangex % column range of position numbers of active wavelets
        rangey % row range of position numbers of active wavelets

        dim % shape of the coefficient matrix, row-by-column
        numel % prod(dim)
        
        pmesh % mesh of integer position index in the range [rangex, rangey]
        pmeshpoints % pmesh reshaped in 2-by-?

        mesh % pmesh * 2^J
        meshpoints % mesh reshaped in 2-by-?
        
        rangex_int % column range of position numbers locating in the interior of the ROI
        rangey_int % row range of position numbers locating in the interior of the ROI

        dim_int 
        numel_int 
        
        pmesh_int % 
        pmeshpoints_int % 

        mesh_int % subset of mesh, the points in the interior of the ROI
        meshpoints_int % mesh_int reshaped in 2-by-?
    end
    
    methods        
        function val = resize(obj, X, rx,ry, mode)
        % val = resize0(obj, X, rx,ry, mode)
        % Resize an image X defined on the grid rx X ry to the current grid according to the
        % specified mode:
        % 'crop': intersection
        % 'xgrid': trunc to the size of X
        % 'ngrid': trunc to the size of current grid
        % 
            val = tools.resizemat(X, rx, ry, obj.rangex, obj.rangey, mode);
        end
        
        function val = resizev(obj, X, V, mode)
        % Resize an image X defined on the grid of the shift-inv space V, to the current
        % space according to the specified mode. Interface function to obj.resize
            
            val = obj.resize(X, V.rangex, V.rangey, mode);
        end
    
        function val = restrict_mask(obj, rx, ry)
        % Make a 0-1 mask for restricting on a smaller mesh defined by rx X ry.
        % INPUTS:
        % rx, ry: index range in x and y axis
        % OUTPUT:
        % val: a 0-1 mask of the same size as the current space, with 1 only on the
        % mesh rx X ry

            if rx(1)<obj.rangex(1) || rx(2)>obj.rangex(2) || ry(1)<obj.rangey(1) || ry(2)> ...
                    obj.rangey(2) % Check dimension first
                error(['Dimension error: the restricted index range must be included in the ' ...
                       'original one.']);
            end
            X = ones(ry(2)-ry(1)+1, rx(2)-rx(1)+1);
            val = tools.resizemat(X, rx, ry, obj.rangex, obj.rangey, 'ngrid');
        end
        
        function val = extend_mask(obj, rx, ry)
        % Make a 0-1 mask for extending to a larger mesh defined by rx X ry.
        % INPUTS:
        % rx, ry: index range in x and y axis
        % OUTPUT:
        % val: a 0-1 mask of the same size as the grid of rx X ry, with 1 only on
        % the current space

            if rx(1)>obj.rangex(1) || rx(2)<obj.rangex(2) || ry(1)>obj.rangey(1) || ry(2)< ...
                    obj.rangey(2) % Check dimension first
                error(['Dimension error: the extended index range must include the original ' ...
                       'one.']);                
            end
            X = ones(ry(2)-ry(1)+1, rx(2)-rx(1)+1);
            val = tools.resizemat(X, rx, ry, obj.rangex, obj.rangey, 'xgrid');
        end

        function val = interior(obj, X)
        % Restrict X on the interior mesh
            val = tools.resizemat(X, obj.rangex, obj.rangey, obj.rangex_int, obj.rangey_int, 'ngrid');
        end
        
        function val = get.pmesh(obj)
            Nx = obj.rangex(1):obj.rangex(2);
            Ny = obj.rangey(1):obj.rangey(2);
            [val.Sx, val.Sy] = meshgrid(Nx, Ny);
        end

        function val = get.pmeshpoints(obj)
            val = [obj.pmesh.Sx(:) obj.pmesh.Sy(:)]';
        end
        
        function val = get.mesh(obj)
            val.Sx = 2^obj.J * obj.pmesh.Sx;
            val.Sy = 2^obj.J * obj.pmesh.Sy;
        end

        function val = get.meshpoints(obj)
        % val = [obj.mesh.Sx(:) obj.mesh.Sy(:)]';
            val = [obj.pmesh.Sx(:) obj.pmesh.Sy(:)]' * 2^obj.J;
        end

        function val = get.rangex_int(obj)
            rangex1 = ceil(obj.suppx(1)/2^obj.J);
            rangex2 = floor(obj.suppx(2)/2^obj.J);
            val = [rangex1 rangex2];
        end
        
        function val = get.rangey_int(obj)
            rangey1 = ceil(obj.suppy(1)/2^obj.J);
            rangey2 = floor(obj.suppy(2)/2^obj.J);
            val = [rangey1 rangey2];
        end
        
        function val = get.pmesh_int(obj)            
            Nx = obj.rangex_int(1):obj.rangex_int(2);
            Ny = obj.rangey_int(1):obj.rangey_int(2);
            [val.Sx, val.Sy] = meshgrid(Nx, Ny); 
        end

        function val = get.pmeshpoints_int(obj)
            val = [obj.pmesh_int.Sx(:) obj.pmesh_int.Sy(:)]';
        end
        
        function val = get.mesh_int(obj)
            val.Sx = 2^obj.J * obj.pmesh_int.Sx;
            val.Sy = 2^obj.J * obj.pmesh_int.Sy;
        end

        function val = get.meshpoints_int(obj)
        % val = [obj.mesh.Sx(:) obj.mesh.Sy(:)]';
            val = [obj.pmesh_int.Sx(:) obj.pmesh_int.Sy(:)]' * 2^obj.J;
        end
        
        
        function val = get.dim(obj)
            ncol = obj.rangex(2) - obj.rangex(1) + 1;
            nrow = obj.rangey(2) - obj.rangey(1) + 1;
            val = [nrow, ncol]; % row-col
        end
        
        function val = get.numel(obj)
            val = prod(obj.dim);
        end
        
        function val = get.dim_int(obj)
            ncol = obj.rangex_int(2) - obj.rangex_int(1) + 1;
            nrow = obj.rangey_int(2) - obj.rangey_int(1) + 1;
            val = [nrow, ncol]; % row-col
        end
        
        function val = get.numel_int(obj)
            val = prod(obj.dim_int);
        end
        
        function obj = ShiftInvSpace(center, width, height, J, id, qmf_range)
        % obj = ShiftInvSpace(center, width, height, J, id, N1, N2)
        % Construct a shift invariant space on a ROI.
        % INPUTS:
        % qmf_range: if non empty, they define the qmf index range and the orthogonal
        % basis is used. For example, [0, 3] for Daubechies 2 wavelet.
            
            obj = obj@wavelet.ROI(center, width, height);
            obj.J = J;
            obj.id = id;
            obj.qmf_range = qmf_range;
            
            if ~isempty(qmf_range)
                obj.ortho = 1;
            else
                obj.ortho = 0;
            end
            
            if obj.ortho
                [obj.rangex, obj.rangey] = wavelet.OrthoWvl.active_position(obj.center, obj.width, ...
                                                                  obj.height, qmf_range(1), ...
                                                                  qmf_range(2), obj.J, obj.id);
            else
                error('Not implemented: ShiftInvSpace for redundant basis.');
            end
        end        
    end
    
end


% function val = extend0(obj, X, nrx, nry)
% % val = extend0(obj, X, nrx, nry)
% %
% % Extend (by zero-padding) a matrix X defined on the current mesh to a larger mesh where the
% % range of positions are nrx and nry.

%     if size(X,1)~=obj.dim(1) || size(X,2)~=obj.dim(2)
%         error(['Dimension error: the input matrix must have same dimension as obj.dim']);
%     end

%     rx = obj.rangex; ry = obj.rangey;
%     if rx(1)<nrx(1) || rx(2)>nrx(2) || ry(1)<nry(1) || ry(2)> ...
%             nry(2) % Check dimension first
%         error(['Dimension error: the extended index range must include the original ' ...
%                'one.']);
%     end
%     val = zeros(nry(2)-nry(1)+1, nrx(2)-nrx(1)+1);
%     col1 = rx(1)-nrx(1)+1; col2 = col1+rx(2)-rx(1);
%     row1 = ry(1)-nry(1)+1; row2 = row1+ry(2)-ry(1);            
%     val(row1:row2, col1:col2) = X;
% end

% function val = restrict0(obj, X, rx, ry)
% % val = restrict0(obj, X, rx, ry)
% %
% % Restrict (by zero-padding) a matrix X defined on the current mesh to a smaller mesh where
% % the range of positions are rx and ry.

%     if size(X,1)~=obj.dim(1) || size(X,2)~=obj.dim(2)
%         error(['Dimension error: the input matrix must have same dimension as obj.dim']);
%     end

%     nrx = obj.rangex; nry = obj.rangey;
%     if rx(1)<nrx(1) || rx(2)>nrx(2) || ry(1)<nry(1) || ry(2)> ...
%             nry(2) % Check dimension first
%         error(['Dimension error: the restricted index range must be included the original ' ...
%                'one.']);
%     end
%     val = zeros(size(X));
%     col1 = rx(1)-nrx(1)+1; col2 = col1+rx(2)-rx(1);
%     row1 = ry(1)-nry(1)+1; row2 = row1+ry(2)-ry(1);            
%     val(row1:row2, col1:col2) = X(row1:row2, col1:col2);
% end

% function val = crop0(obj, X, rx, ry)
% % val = crop0(obj, X, rx, ry)
% % Crop a matrix X defined on the current mesh to another mesh where the
% % range of positions are rx and ry.

%     row1 = max(ry(1), obj.rangey(1));
%     row2 = min(ry(2), obj.rangey(2));
%     col1 = max(rx(1), obj.rangex(1));
%     col2 = min(rx(2), obj.rangex(2));

%     mask = obj.restrict0(ones(obj.dim), [col1, col2], [row1, row2]);
%     val = reshape(X(find(mask)), row2-row1+1, col2-col1+1);
% end

% function val = extend(obj, X, V)
% % val = extend(obj, X, V)
% % Interface to the function extend0, with V another ShiftInvSpace instance
%     val = obj.extend0(X, V.rangex, V.rangey);
% end

% function val = restrict(obj, X, V)
% % Interface to the function restrict0, with V another ShiftInvSpace instance
%     val = obj.restrict0(X, V.rangex, V.rangey);
% end

% function val = crop(obj, X, V)
% % Interface to the function crop0, with V another ShiftInvSpace instance
%     val = obj.crop0(X, V.rangex, V.rangey);
% end
