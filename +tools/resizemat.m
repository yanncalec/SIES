function val = resizemat(X, rx, ry, nrx, nry, mode)
% Resize an image X defined on the grid rx X ry to the another grid nrx X nry,
% according to the specified mode.
% 'crop': intersection
% 'xgrid': trunc to the grid rx X ry
% 'ngrid': trunc to the grid nrx X nry
% 

    [nrow, ncol] = size(X);
    if nrow ~= ry(2)-ry(1)+1 || ncol ~= rx(2)-rx(1)+1
        error('Size of X and the grid mismatch');
    end
    
    ry1 = max(ry(1), nry(1));
    ry2 = min(ry(2), nry(2));
    rx1 = max(rx(1), nrx(1));
    rx2 = min(rx(2), nrx(2));
    
    dim = [ry2-ry1+1, rx2-rx1+1]; % dimension of the intersection

    ndim = [nry(2)-nry(1)+1, nrx(2)-nrx(1)+1]; % dimension of the new grid
    
    % Starting index of the intersection in the grid nrx X nry
    row_g = ry1-nry(1)+1; 
    col_g = rx1-nrx(1)+1;

    % Starting index of the intersection in the grid of rx X ry
    row_x = ry1-ry(1)+1; 
    col_x = rx1-rx(1)+1;

    if strcmp(mode, 'crop')
        val = X(row_x:row_x+dim(1)-1, col_x:col_x+dim(2)-1); 
    elseif strcmp(mode, 'ngrid')
        val = zeros(ndim);
        val(row_g:row_g+dim(1)-1, col_g:col_g+dim(2)-1) = X(row_x:row_x+dim(1)-1, ...
                                                          col_x:col_x+dim(2)-1); 
    elseif strcmp(mode, 'xgrid')
        val = zeros(size(X));
        val(row_x:row_x+dim(1)-1, col_x:col_x+dim(2)-1) = X(row_x:row_x+dim(1)-1, ...
                                                          col_x:col_x+dim(2)-1); 
    end
end
