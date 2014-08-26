function [D,t] = boundarydet(acc, IM)
% [D,t] = boundarydet(acc, IM, Np)
%
% This function is borrowed from Y.Capdeboscq and A.B.Karmann's
% code. Original help message:
%
% "This function parameterizes a simply connected region from an image.
% Make the image in paint or some other imaging program.  Black
% foreground and white background works well.  The acc parameter is the
% roughly the distance between each consecutive pair of points in the
% parameterization (in image pixels).  IM is a matrix of the image
% once the image has been discretized into binary values (and thus 
% separated into a foreground and background).  We attribute the level
% set routines to Karrman and Allaire (2009).  A quick description of
% the algorithm is that we start at a point on the boundary, advect a
% point along the tangent for a distance equal to acc, then return to 
% the boundary using the normal vector (defined from the level-set).
% We repeat this process until until we return to the original point in
% the parameterization.
%  
%  Created 23-SEP-2011 by  Anton Bongio Karrman karrman@caltech.edu"
%
% We have simply modified the oridinal function such that the final
% boundary D and t are not tied off.

% We measure the space the image is in by pixels:
    [ny,nx] = size(IM) ;
    dx = 1 ;
    dy = dx ;

    % We use a level-set routine to define a level-set for the image
    % boundary (This level-set is key in the parameterization algorithm):
    
    % We first initialize the level set with RIiter number of iterations:
    RIiter = 20 ;   
    [X,Y] = meshgrid(1:nx,1:ny) ;  
    
    phi = 2*(IM-.5) ;      
    cfl = 0.5 ;   % CFL CONDITION
                  % Using our CFL condition, we define the time step as
    dt0 = min(dx,dy)*cfl ;

    for n = 1 : RIiter                

        % For our first derivatives:
        phin = shift2n('n',phi) ;
        phis = shift2n('s',phi) ;
        phie = shift2n('e',phi) ;
        phiw = shift2n('w',phi) ;

        % Our scheme is second order in space, so
        % we use these for our second derivatives:
        phinn = shift2n('n',phin) ;
        phiss = shift2n('s',phis) ;
        phiee = shift2n('e',phie) ;
        phiww = shift2n('w',phiw) ;

        %  The first derivatives:
        dxm = (phi-phiw)/dx ;
        dxp = (phie-phi)/dx ;
        dym = (phi-phis)/dy ;
        dyp = (phin-phi)/dy ;

        % The second derivatives (our scheme is second
        % order so we use several different ones):
        dxmxm = (phi - 2*phiw + phiww)/(dx^2) ;
        dxpxp = (phiee - 2*phie + phi)/(dx^2) ;
        dxpxm = (phie - 2*phi + phiw)/(dx^2) ;

        dymym = (phi - 2*phis + phiss)/(dy^2) ;
        dypyp = (phinn - 2*phin + phi)/(dy^2) ;
        dypym = (phin - 2*phi + phis)/(dy^2) ;

        % From Sethian, our scheme requires four
        % main parts, as defined here:
        partA = dxm + .5*dx*minmod(dxmxm,dxpxm) ;
        partB = dxp - .5*dx*minmod(dxpxp,dxpxm) ;
        partC = dym + .5*dy*minmod(dymym,dypym) ;
        partD = dyp - .5*dy*minmod(dypyp,dypym) ;

        % Then we use these parts along with the 
        % flux function to get:
        delp2 = g(partA,partB,partC,partD) ;
        delm2 = g(partB,partA,partD,partC) ;

        % sphi is the sign of phi and we added a small 
        % value proportional to the mesh size in the 
        % denominator to be sure we never divide by zero.
        nabla = 0.5*(dxm.^2 + dxp.^2 + dym.^2 + dyp.^2) ;
        sphi = phi./sqrt(phi.*phi+sqrt(dx^2+dy^2)*nabla/10) ;
        sphip = max(sphi,0) ;
        sphim = min(sphi,0) ;

        phi = phi - dt0*(sphip.*delp2 + sphim.*delm2 - sphi) ;
        %evalin('caller','set(handles.paramstatus,''String'',num2str(round(toc)))') ;   
    end
    % Once our iterations are finished, we set 
    % our output to the new level set function:
    phix = (phie-phiw)/(2*dx) ;
    phiy = (phin-phis)/(2*dy) ;
    % We define an epsilon value for the level-set definitions:
    epscurv = min(dx,dy)/40 ;
    % We define some level-set terms:
    % Norm of gradient:
    mag = sqrt(phix.^2+phiy.^2+epscurv^2) ;    
    % X, Y component of normal, tangent vector
    normx = phix./mag ; 
    normy = phiy./mag ;   
    tanx = -normy ;
    tany = normx ;    
    
% $$$     hold on ;  
% $$$     contour(X,Y,IM) ;

    % The beginning of the algorithm starts with 2 runs to make sure we
    % start as close to the boundary as possible:  
    
    % We start at the point in the shape halfway up the image, farthest
    % to the right:
    ny2 = interp1(Y(:,1),Y(:,1),ny/2,'nearest') ;        
    x0 = [max(find(round(phi(Y == ny2))<0)) round(ny2)] ;    
    
    for i = 1:2
        % We evaluate the tangent vector as the x0:
        txy = [interp2(X,Y,tanx,x0(1),x0(2)) interp2(X,Y,tany,x0(1),x0(2))] ;
        % We advect our point, a distance equal to the acc parameter, along
        % the tangent vector:
        x0 = x0+acc*txy ;
        % Once we have shot our point off the tangent direction, we then
        % return to the boundary by following the normal vector (the
        % distance we travel is the distance away from the
        % boundary, which is the signed distance function, which is the
        % same thing as the level-set function phi):
        nxy = [interp2(X,Y,normx,x0(1),x0(2)) interp2(X,Y,normy,x0(1),x0(2))] ;   
        x0 = x0-interp2(X,Y,phi,x0(1),x0(2))*nxy ;   
    end        
    
    % Once we are at the boundary, we enter into the loop that finds the
    % entire parameterization:    
    D = x0 ;
    k = 1 ;
    % After we have gotten a few steps (4 to be precise) into the
    % algorithm, we end the parameterization when we have gone around the
    % perimeter and gotten close (less than 2*acc "pixels") to the original
    % x0:
    while (norm(D(1,:)-D(end,:)) > 2*acc || k<4)
        % We start at the last point in the parameterization:
        x = D(end,:) ;
        % We advect along the tangent a distance equal to the acc
        % parameter:
        txy = [interp2(X,Y,tanx,x(1),x(2)) interp2(X,Y,tany,x(1),x(2))] ;
        x = x+acc*txy ;  
        % We return to the boundary by advecting along the normal vector a
        % distance given by the signed-distance function (phi, the level
        % set):
        nxy = [interp2(X,Y,normx,x(1),x(2)) interp2(X,Y,normy,x(1),x(2))] ;
        x = x-interp2(X,Y,phi,x(1),x(2))*nxy ;                     
        D(end+1,:) = x ;      
        k = k+1 ;
% $$$         % Print the time the algorithm is taking in the parameterization
% $$$         % status window:
% $$$         evalin('caller','set(handles.paramstatus,''String'',strcat(''Time Passed:'',num2str(round(toc))))') ;      
% $$$         drawnow ;        
% $$$         % If we click the cancel button, we cancel the parameterization: 
% $$$         if get(gca,'UserData') == 1
% $$$             break ;
% $$$         end              
    end
    
    %      % We redefine the parameter and tie off the the paramerization:
    %      t = (0:2*pi/length(D):2*pi) ;      % length of t is length(D)+1
    %      D = [D;D(1,:)]' ;

    % By definition, D is a closed curve. But we remove the end
    % point of D and t such that the first and the last elements in
    % D and t are not the same (not tied-off). 
    t = 2*pi*(0:length(D)-1)/length(D) ;
    D = D' ; 
    
    %     % We interpolate D wtih 1024 (2^10) points (this is a fine
    %     % interpolation that can be reduced later).
    %     tnew = 2*pi*(0:Np-1)/Np;
    %     D = interp1(t,D,tnew,'spline') ;
    %     D = D';
    %     t = tnew;
end

% SPACE SHIFT FUNCTION
% Using Neumann boundary conditions, we just shift
% our matrix over one index in a certain direction
% in order to take derivatives. 
function phishift = shift2n(direction,phi)
% SHIFTS LEVEL SET FUNCTION WITH NEUMANN CONDITIONS
    switch direction
      case 'e'           % SHIFT WEST
        [m,n] = size(phi) ;
        phishift(1:m,1:n-1) = phi(1:m,2:n) ;
        phishift(1:m,n) = phi(1:m,n) ;
      case 'w'           % SHIFT EAST
        [m,n] = size(phi) ;
        phishift(1:m,2:n) = phi(1:m,1:n-1) ;
        phishift(1:m,1) = phi(1:m,1) ;
      case 'n'           % SHIFT NORTH
        [m,n] = size(phi) ;
        phishift(1:m-1,1:n) = phi(2:m,1:n) ;
        phishift(m,1:n) = phi(m,1:n) ;
      case 's'           % SHIFT SOUTH
        [m,n] = size(phi) ;
        phishift(2:m,1:n) = phi(1:m-1,1:n) ;
        phishift(1,1:n) = phi(1,1:n) ;
    end
end        

% MINMOD FUNCTION
% This function is essential for our second order
% scheme used to solve the level set equation 
% (as well as for for reinitialization).
function mout = minmod(phi1,phi2)
    sphi1=sign(phi1) ;
    sphi2=sign(phi2) ;
    mout = max(0,sphi1.*sphi2).*sphi1.*min(abs(phi1),abs(phi2)) ;
end

% FLUX FUNCTION

% This is the numerical flux of our scheme for Hamilton Jacobi equations
function gout = g(u1,u2,v1,v2)
    gout = sqrt( max(u1,0.).^2 + min(u2,0.).^2 ...
                 + max(v1,0.).^2 + min(v2,0.).^2 ) ;
end

%  
%  Created 23-SEP-2011 by  Anton Bongio Karrman karrman@caltech.edu
%  Further credit to Gregoire Allaire for the level-set routines
%

