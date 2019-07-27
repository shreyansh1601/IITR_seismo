% Function Poisson_smoother()
% This function makes specified number (iternum) of Gauss Seidel iteraions 
% using relaxation koefficients (krelaxs,krelaxc)
% for Stokes and Continuity equations defined on 2D staggered grid
% with specified resolution (xnum, ynum) and gridsteps (xstp, ystp)
% given distribution of right parts for all equations (RX,RY,RC) on the grid 
% and given viscosity (etas) on the grid 
% pressure is normalized relative to given value (prnorm) in the first cell 
%
% Function return new approximation for velocity and pressure (vx,vy,pr)
% and distribution of residuals (resx,resy,resc)
function[vx,resx,vy,resy,pr,resc]=Stokes_Continuity_smoother(prnorm,etas,iternum,krelaxs,krelaxc,xnum,ynum,xstp,ystp,RX,vx,RY,vy,RC,pr)
% 
% Staggered Grid for Multigrid
% 
%     vx       vx       vx    
%
% vy  +---vy---+---vy---+   vy
%     |        |        |
%     vx   P   vx   P   vx    
%     |        |        |
% vy  +---vy---+---vy---+   vy
%     |        |        |
%     vx   P   vx   P   vx    
%     |        |        |
% vy  +---vy---+---vy---+   vy
%
%     vx       vx       vx    
% 
% Lines show basic grid
% Basic (density) nodes are shown with +
% Ghost nodes shown outside the basic grid
% are used for boundary conditions


% Poisson-like equations koefficients
xkf=etas/xstp^2;
ykf=etas/ystp^2;
xykf=-etas*(2/xstp^2+2/ystp^2);

% Gauss-Seidel vx, xy, P iteration circle
for niter=1:1:iternum;
    % Solving of Stokes and continuity equations on nodes
    for i=1:1:ynum+1;
        for j=1:1:xnum+1;

            % x-Stokes equation d2vx/dx2+d2vx/dy2-dP/dx=RX
            if (j<xnum+1)
                % vx-Boundrary conditions 
                if (i==1 | i==ynum+1 | j==1 | j==xnum);
                   % Boundary condition for vx
                    vx(i,j)=0;
                    % Upper Boundary
                    if(i==1)
                        % Free slip
%                        vx(i,j)=vx(i+1,j);
                        % No slip
                        vx(i,j)=-vx(i+1,j);
                    end
                    % Lower boundary
                    if(i==ynum+1)
                        % Free slip
%                        vx(i,j)=vx(i-1,j);
                        % No slip
                        vx(i,j)=-vx(i-1,j);
                    end
                % Solving x-Stokes equation
                else
                    % Computing current residual
                    resx(i,j)=RX(i,j)-(xkf*(vx(i,j-1)+vx(i,j+1))+ykf*(vx(i-1,j)+vx(i+1,j))+xykf*vx(i,j)-(pr(i-1,j)-pr(i-1,j-1))/xstp);
                    % Updating solution
                    vx(i,j)=vx(i,j)+resx(i,j)/xykf*krelaxs;
                end
            end
            
            % y-Stokes equation d2vy/dx2+d2vy/dy2-dP/dy=RY
            if (i<ynum+1)
                % vy-Boundrary conditions 
                if (i==1 | i==ynum | j==1 | j==xnum+1);
                   % Boundary condition for vx
                    vy(i,j)=0;
                    % Left boundary
                    if(j==1)
                        % Free slip
                        vy(i,j)=vy(i,j+1);
                        % No slip
%                        vy(i,j)=-vy(i,j+1);
                    end
                    % Right boundary
                    if(j==xnum+1)
                        % Free slip
                        vy(i,j)=vy(i,j-1);
                        % No slip
%                        vy(i,j)=-vy(i,j-1);
                    end
                % Solving y-Stokes equation
                else
                    % Computing current residual
                    resy(i,j)=RY(i,j)-(xkf*(vy(i,j-1)+vy(i,j+1))+ykf*(vy(i-1,j)+vy(i+1,j))+xykf*vy(i,j)-(pr(i,j-1)-pr(i-1,j-1))/ystp);
                    % Updating solution
                    vy(i,j)=vy(i,j)+resy(i,j)/xykf*krelaxs;
                end
            end
            
            % Continuity equation dvx/dx+dvy/dy=RC
            % is solved via pressure updates
            % dpr=-etas*div(v)
            % pr-Boundrary conditions 
            if (i<ynum && j<xnum);
            % Solving Continuity equation by adjusting pressure
                % Computing current residual
                resc(i,j)=RC(i,j)-((vx(i+1,j+1)-vx(i+1,j))/xstp+(vy(i+1,j+1)-vy(i,j+1))/ystp);
                % Updating pressure solution
                pr(i,j)=pr(i,j)+resc(i,j)*etas*krelaxc;
            end
             
        end            
    end
    % End of Solving Stokes and continuity equations on nodes
end

% Correct pressure with the given value in the first cell
dp=prnorm-pr(1,1);
pr=pr+dp;

% Computing final state of residuals
for i=1:1:ynum+1;
    for j=1:1:xnum+1;

        % x-Stokes equation d2vx/dx2+d2vx/dy2-dP/dx=RX
        if (j<xnum+1)
            % vx-Boundrary conditions 
            if (i==1 | i==ynum+1 | j==1 | j==xnum);
                resx(i,j)=0;
            % x-Stokes equation
            else
                % Computing current residual
                resx(i,j)=RX(i,j)-(xkf*(vx(i,j-1)+vx(i,j+1))+ykf*(vx(i-1,j)+vx(i+1,j))+xykf*vx(i,j)-(pr(i-1,j)-pr(i-1,j-1))/xstp);
            end
        end
            
        % y-Stokes equation d2vy/dx2+d2vy/dy2-dP/dy=RY
        if (i<ynum+1)
            % vy-Boundrary conditions 
            if (i==1 | i==ynum | j==1 | j==xnum+1);
                resy(i,j)=0;
            %y-Stokes equation
            else
                % Computing current residual
                resy(i,j)=RY(i,j)-(xkf*(vy(i,j-1)+vy(i,j+1))+ykf*(vy(i-1,j)+vy(i+1,j))+xykf*vy(i,j)-(pr(i,j-1)-pr(i-1,j-1))/ystp);
            end
        end
            
        % Continuity equation dvx/dx+dvy/dy=RC
        % is solved via pressure updates
        % dpr=-etas*div(v)
        if (i<ynum && j<xnum);
            % Computing current residual
            resc(i,j)=RC(i,j)-((vx(i+1,j+1)-vx(i+1,j))/xstp+(vy(i+1,j+1)-vy(i,j+1))/ystp);
        end
             
    end            
end


