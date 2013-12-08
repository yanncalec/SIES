function plot_field(obj, F, F_bg, SX, SY, nbLine, varargin)

% First, troncate the infinite terms
    matrix_isinf = isinf(F) ;
    F(matrix_isinf) = mean(mean(F(isfinite(F)))) ;

    matrix_isinf = isinf(F_bg) ;
    F_bg(matrix_isinf) = mean(mean(F_bg(isfinite(F_bg)))) ;

    % Field could be complex, plot the real part
    figure; 
    subplot(2,2,1);
    contourf(SX, SY, real(F), nbLine); hold on;
    for n=1:obj.nbIncls        
        plot(obj.D{n}, varargin{:});
    end
    axis image; 
    title('Field u, real part');
    colorbar();

    subplot(2,2,2);
    contourf(SX, SY, imag(F), nbLine); hold on;
    for n=1:obj.nbIncls        
        plot(obj.D{n}, varargin{:});
    end
    axis image; 
    title('Field u, imaginary part');
    colorbar();

    subplot(2,2,3);
    contourf(SX, SY, real(F-F_bg), nbLine); hold on;
    for n=1:obj.nbIncls        
        plot(obj.D{n}, varargin{:});
    end
    axis image;
    title('Perturbed field u-U, real part');
    colorbar();

    % Background field is real
    subplot(2,2,4);
    contourf(SX, SY, real(F_bg), nbLine); hold on;
    for n=1:obj.nbIncls        
        plot(obj.D{n}, varargin{:});
    end
    axis image; 
    title('Background field U, real part');
    colorbar();

end
