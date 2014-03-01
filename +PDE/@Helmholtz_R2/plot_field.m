function plot_field(obj, s, F, F_bg, SX, SY, nbLine, varargin)

    src = obj.cfg.src(s);
    % Field could be complex, plot the real part

    figure; 
    subplot(2,2,1);
    contourf(SX, SY, real(F), nbLine); hold on;
    plot(src(1), src(2), 'gx', varargin{:});
    for i=1:obj.nbIncls
        plot(obj.D{i}, varargin{:});
    end
    axis image; 
    title('Field u, real part');
    colorbar();

    subplot(2,2,2);
    contourf(SX, SY, imag(F), nbLine); hold on;
    plot(src(1), src(2), 'gx', varargin{:});
    for i=1:obj.nbIncls
        plot(obj.D{i}, varargin{:});
    end
    axis image; 
    title('Field u, imaginary part');
    colorbar();

    subplot(2,2,3);
    contourf(SX, SY, real(F-F_bg), nbLine); hold on;
    plot(src(1), src(2), 'gx', varargin{:});
    for i=1:obj.nbIncls
        plot(obj.D{i}, varargin{:});
    end
    axis image;
    title('Perturbed field u-U, real part');
    colorbar();

    % Background field is real
    subplot(2,2,4);
    contourf(SX, SY, real(F_bg), nbLine); hold on;
    plot(src(1), src(2), 'gx', varargin{:});
    % plot(obj.D{1}, varargin{:});
    axis image; 
    title('Background field U, real part');
    colorbar();
    
end