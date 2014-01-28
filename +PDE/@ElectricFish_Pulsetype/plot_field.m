function plot_field(obj, s, F, F_bg, SX, SY, nbLine, subfig, varargin)

% First, troncate the infinite terms
    % matrix_isinf = isinf(F) ;
    % F(matrix_isinf) = mean(mean(F(isfinite(F)))) ;

    % matrix_isinf = isinf(F_bg) ;
    % F_bg(matrix_isinf) = mean(mean(F_bg(isfinite(F_bg)))) ;

    % Field could be complex, plot the real part
    
    figure; 
    
    if subfig
        subplot(2,2,1);
    end
    contourf(SX, SY, real(F), nbLine); hold on;
    plot(obj.cfg.Bodies(s), varargin{:});
    for n=1:obj.nbIncls        
        plot(obj.D{n}, varargin{:});
    end
    axis image; 
    title('Potential field u, real part');
    colorbar();

    if subfig
        subplot(2,2,2);
    else
        figure;
    end
    contourf(SX, SY, imag(F), nbLine); hold on;
    plot(obj.cfg.Bodies(s), varargin{:});
    for n=1:obj.nbIncls        
        plot(obj.D{n}, varargin{:});
    end
    axis image; 
    title('Potential field u, imaginary part');
    colorbar();

    if subfig
        subplot(2,2,3);
    else
        figure;
    end
    contourf(SX, SY, real(F-F_bg), nbLine); hold on;
    plot(obj.cfg.Bodies(s), varargin{:});
    for n=1:obj.nbIncls        
        plot(obj.D{n}, varargin{:});
    end
    axis image;
    title('Perturbed field u-U, real part');
    colorbar();

    % Background field is real
    if subfig
        subplot(2,2,4);
    else
        figure;
    end
    contourf(SX, SY, F_bg, nbLine); hold on;
    plot(obj.cfg.Bodies(s), varargin{:});
    axis image; 
    title('Background potential field U');
    colorbar();    
end
