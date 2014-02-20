function A = Green2D_apprxcoeff(obj, Z)

    Apprx = obj.ApprxSpace{1};

    for nx=Apprx.rangex(1):Apprx.rangex(2) % iteration on column
        for ny=Apprx.rangey(1):Apprx.rangey(2) % iteration on row
            posx = 2^obj.Jmax*nx; posy = 2^obj.Jmax*ny; 
            [rx, ry] = obj.getsupport(obj.Jmax, 0, [nx, ny]');
            nrow = (ny-Apprx.rangey(1)) + 1;
            ncol = nx-Apprx.rangex(1) + 1;

            % Check if the support of Phi contains the singularity of Green function
            if (rx(1)<=Z(1)) && (rx(2)>=Z(1)) && (ry(1)<=Z(2)) && (ry(2)>=Z(2))
                error(['The scaling function of j=',num2str(obj.Jmax), ', n=[', num2str(nx),', ', num2str(ny), '] contains ' ...
                                    'the singularity, try to reduce Jmax.']);
            end
            
            %         SX = tools.tensorplus_mex((obj.Gphi-Z(2)).^2, (obj.Gphi-Z(1)).^2);
            %         Y = 1/4/pi * log(SX);
            %         T = tools.tensorprod_mex(obj.Tphi, obj.Tphi);
            %         A.coeff(nrow, ncol) = sum(sum(Y.*T))*obj.Gphi(1)^2;
        end
    end

    A.coeff =  wavelet.OrthoWvl.Green2D_apprxcoeff_mex(obj.Tphi, obj.phisupp, obj.Jmax, Apprx.rangex, Apprx.rangey, Z);

end

    % X = tools.Laplacian.Green2D(obj.AROI.meshpoints, Z(:)) * 2^obj.Jnum;
    % % decomposition from Jnum to Jmax
    % [~, ~, W] = obj.analysis(reshape(X, obj.AROI.dim)); 
    % % W = W(1:numel(obj.wdim.A));
    % Wa = W(find(obj.wmask.A(:)));
