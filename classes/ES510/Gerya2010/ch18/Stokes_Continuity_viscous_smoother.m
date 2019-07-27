% Function Stokes_Continuity_smoother_ghost()
% This function makes specified number (iternum) of Gauss Seidel iteraions 
% using relaxation koefficients (krelaxs,krelaxc)
% for Stokes and Continuity equations defined on 2D staggered grid
% with specified resolution (xnum, ynum) and gridsteps (xstp, ystp)
% given distribution of right parts for all equations (RX,RY,RC) on the grid 
% and given variable shear (etas) and normal (etan) viscosity distributions 
% pressure is normalized relative to given value (prnorm) in the first cell
%
% Velocity Boundary condition are implemented from ghost nodes 
% directly into Stokes and continuity equations
%
% Function return new approximation for velocity and pressure (vx,vy,pr)
% and distribution of residuals (resx,resy,resc)
function[vx,resx,vy,resy,pr,resc]=Stokes_Continuity_viscous_smoother(prnorm,etas,etan,iternum,krelaxs,krelaxc,xnum,ynum,xstp,ystp,RX,vx,RY,vy,RC,pr)
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
xkf=1/xstp^2;
xkf2=2/xstp^2;
ykf=1/ystp^2;
ykf2=2/ystp^2;
xykf=1/(xstp*ystp);

% Gauss-Seidel vx, xy, P iteration circle
for niter=1:1:iternum;
    % Solving of Stokes and continuity equations on nodes
    for i=1:1:ynum+1;
        for j=1:1:xnum+1;

            % x-Stokes equation d2vx/dx2+d2vx/dy2-dP/dx=RX
            if (j<xnum+1)
                % Solving x-Stokes equation 
                if (i>1 && i<ynum+1 && j>1 && j<xnum);
                    
                    % Initial value of koefficient for normal internal nodes
                    kfxcur=-xkf2*(etan(i-1,j)+etan(i-1,j-1))-ykf*(etas(i,j)+etas(i-1,j));
                    % Marginal nodes
                    % Left boundary nodes
                    if (j==2)
                        %No slip, Free slip
                        vx(i,j-1)=0;
                        %kfxcur=kfxcur-0;
                    else
                        % Right boundary nodes
                        if (j==xnum-1)
                            %No slip, Free slip
                            vx(i,j+1)=0;
                            %kfxcur=kfxcur+0;
                        end
                    end
                    % Upper boundary nodes
                    if (i==2)
                        %Free slip
                        vx(i-1,j)=vx(i,j);
                        kfxcur=kfxcur+etas(i-1,j)*ykf;
                        %No slip
%                       vx(i-1,j)=-vx(i,j);
%                       kfxcur=kfxcur-etas(i-1,j)*ykf;
                    else
                        % Lower boundary nodes
                        if (i==ynum)
                            %Free slip
                            vx(i+1,j)=vx(i,j);
                            kfxcur=kfxcur+etas(i,j)*ykf;
                            %No slip: vx(i,j+1)=-vx(i,j)
%                           vx(i+1,j)=-vx(i,j);
%                           kfxcur=kfxcur-etas(i,j)*ykf;
                        end
                    end
                    % x-Stokes equation stensil
                    %     +-------------------- -+----------------------+   
                    %     |                      |                      |
                    %     |                      |                      |
                    %     |                   vx(i-1,j)                 |    
                    %     |                      |                      |
                    %     |                      |                      |
                    %     +-----vy(i-1,j)---etas(i-1,j)---vy(i-1,j+1)---+
                    %     |                      |                      |
                    %     |                      |                      |
                    % vx(i,j-1)  pr(i-1,j-1)  vx(i,j)     P(i-1,j)   vx(i,j+1)    
                    %     |     etan(i-1,j-1)    |       etan(i-1,j)    |
                    %     |                      |                      |
                    %     +------vy(i,j)-----etas(i,j)----vy(i,j+1)-----+
                    %     |                      |                      |
                    %     |                      |                      |
                    %     |                   vx(i+1,j)                 |    
                    %     |                      |                      |
                    %     |                      |                      |
                    %     +-------------------- -+----------------------+   
                    % Computing Current x-Stokes residual
                    % dSIGMAxx/dx-dP/dx
                    resxcur=RX(i,j)-(xkf2*(etan(i-1,j)*(vx(i,j+1)-vx(i,j))-etan(i-1,j-1)*(vx(i,j)-vx(i,j-1)))-(pr(i-1,j)-pr(i-1,j-1))/xstp);
                    % dSIGMAxy/dy
                    resxcur=resxcur-(etas(i,j)*((vx(i+1,j)-vx(i,j))*ykf+(vy(i,j+1)-vy(i,j))*xykf)-etas(i-1,j)*((vx(i,j)-vx(i-1,j))*ykf+(vy(i-1,j+1)-vy(i-1,j))*xykf));
                    % Updating solution
                    vx(i,j)=vx(i,j)+resxcur/kfxcur*krelaxs;

                end
            end
            
            
            % y-Stokes equation d2vy/dx2+d2vy/dy2-dP/dy=RY
            if (i<ynum+1)
                % Solving y-Stokes equation
                if (i>1 && i<ynum && j>1 && j<xnum+1);
                    % Initial value of koefficient for normal internal nodes
                    kfycur=-ykf2*(etan(i,j-1)+etan(i-1,j-1))-xkf*(etas(i,j)+etas(i,j-1));
                    % Marginal nodes
                    % Left boundary nodes
                    if (j==2)
                        %Free slip
                        vy(i,j-1)=vy(i,j);
                        kfycur=kfycur+etas(i,j-1)*xkf;
                        %No slip
%                       vy(i,j-1)=-vy(i,j);
%                       kfycur=kfycur-etas(i,j-1)*xkf;
                    else
                        % Right boundary nodes
                        if (j==xnum)
                            %Free slip
                            vy(i,j+1)=vy(i,j);
                            kfycur=kfycur+etas(i,j)*xkf;
                            %No slip
%                           vy(i,j+1)=-vy(i,j);
%                           kfycur=kfycur-etas(i,j)*xkf;
                        end
                    end
                    % Upper boundary nodes
                    if (i==2)
                        %No sleep Free slip 
                        vy(i-1,j)=0;
                        %kfycur=kfycur+0;
                    else
                        % Lower boundary nodes
                        if (i==ynum-1)
                            %No sleep Free slip
                            vy(i+1,j)=0;
                            %kfxcur=kfxcur+0;
                        end
                    end
                    % y-Stokes equation stensil
                    %     +-------------------- -+-------vy(i-1,j)------+----------------------+    
                    %     |                      |                      |                      |
                    %     |                      |                      |                      |
                    %     |                  vx(i,j-1)   P(i-1,j-1)  vx(i,j)                   |    
                    %     |                      |      etan(i-1,j-1)   |                      |
                    %     |                      |                      |                      |
                    %     +-----vy(i,j-1)---etas(i,j-1)---vy(i,j)--etas(i,j)-----vy(i,j+1)-----+
                    %     |                      |                      |                      |
                    %     |                      |                      |                      |
                    %     |                  vx(i+1,j-1)  P(i,j-1)   vx(i+1,j)               |    
                    %     |                      |      etan(i,j-1)     |                      |
                    %     |                      |                      |                      |
                    %     +----------------------+-------vy(i+1,j)------+----------------------+
                    %
                    % Computing Current y-Stokes residual
                    % dSIGMAyy/dy-dP/dy
                    resycur=RY(i,j)-(ykf2*(etan(i,j-1)*(vy(i+1,j)-vy(i,j))-etan(i-1,j-1)*(vy(i,j)-vy(i-1,j)))-(pr(i,j-1)-pr(i-1,j-1))/ystp);
                    % dSIGMAxy/dx
                    resycur=resycur-(etas(i,j)*((vy(i,j+1)-vy(i,j))*xkf+(vx(i+1,j)-vx(i,j))*xykf)-etas(i,j-1)*((vy(i,j)-vy(i,j-1))*xkf+(vx(i+1,j-1)-vx(i,j-1))*xykf));
                    % Updating solution
                    vy(i,j)=vy(i,j)+resycur/kfycur*krelaxs;
                    
                end
            end
            
            % Continuity equation dvx/dx+dvy/dy=RC
            % is solved via pressure updates
            % dpr=-etas*div(v)
            % pr-Boundrary conditions 
            if (i<ynum && j<xnum);
            % Solving Continuity equation by adjusting pressure
                % Marginal nodes
                % Left boundary nodes
                if (j==1)
                    %Free slip, No slip
                    vx(i+1,j)=0;
                else
                    % Right boundary nodes
                    if (j==xnum-1)
                        %Free slip, No slip
                        vx(i+1,j+1)=0;
                    end
                end
                % Upper boundary nodes
                if (i==1)
                    %No sleep, Free slip 
                    vy(i,j+1)=0;
                else
                    % Lower boundary nodes
                    if (i==ynum-1)
                        %No sleep, Free slip 
                        vy(i+1,j+1)=0;
                    end
                end
                % Computing current residual
                resccur=RC(i,j)-((vx(i+1,j+1)-vx(i+1,j))/xstp+(vy(i+1,j+1)-vy(i,j+1))/ystp);
                % Updating pressure solution by using normal viscosity
                pr(i,j)=pr(i,j)+resccur*etan(i,j)*krelaxc;
            end
             
        end            
    end
    % End of Solving Stokes and continuity equations on nodes
end

% Correct pressure with the given value in the first cell
dp=prnorm-pr(1,1);
pr=pr+dp;

% Apply vx boundary conditions
% Upper Boundary
% Free slip
 vx(1,:)=vx(2,:);
% No slip
%vx(1,:)=-vx(2,:);
% Lower Boundary
% Free slip
 vx(ynum+1,:)=vx(ynum,:);
% No slip
%vx(ynum+1,:)=-vx(ynum,:);
% Left Boundary
% No slip, Free slip
vx(:,1)=0;
% Right Boundary
% No slip, Free slip
vx(:,xnum)=0;

% Apply vy boundary conditions
% Upper Boundary
% No Sleep, Free slip
 vy(1,:)=0;
% Lower Boundary
% No Sleep, Free slip
 vy(ynum,:)=0;
% Left Boundary
% Free slip
vy(:,1)=vy(:,2);
% No slip
%vy(:,1)=-vy(:,2);
% Right Boundary
% Free slip
vy(:,xnum+1)=vy(:,xnum);
% No slip
%vy(:,xnum+1)=-vy(:,xnum);



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
                % Computing Current x-Stokes residual
                % dSIGMAxx/dx-dP/dx
                resx(i,j)=RX(i,j)-(xkf2*(etan(i-1,j)*(vx(i,j+1)-vx(i,j))-etan(i-1,j-1)*(vx(i,j)-vx(i,j-1)))-(pr(i-1,j)-pr(i-1,j-1))/xstp);
                % dSIGMAxy/dy
                resx(i,j)=resx(i,j)-(etas(i,j)*((vx(i+1,j)-vx(i,j))*ykf+(vy(i,j+1)-vy(i,j))*xykf)-etas(i-1,j)*((vx(i,j)-vx(i-1,j))*ykf+(vy(i-1,j+1)-vy(i-1,j))*xykf));
%                resx(i,j)=RX(i,j)-(xkf1*(vx(i,j-1)+vx(i,j+1))+ykf1*(vx(i-1,j)+vx(i+1,j))+xykf1*vx(i,j)-(pr(i-1,j)-pr(i-1,j-1))/xstp);
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
                % dSIGMAyy/dy-dP/dy
                resy(i,j)=RY(i,j)-(ykf2*(etan(i,j-1)*(vy(i+1,j)-vy(i,j))-etan(i-1,j-1)*(vy(i,j)-vy(i-1,j)))-(pr(i,j-1)-pr(i-1,j-1))/ystp);
                % dSIGMAxy/dx
                resy(i,j)=resy(i,j)-(etas(i,j)*((vy(i,j+1)-vy(i,j))*xkf+(vx(i+1,j)-vx(i,j))*xykf)-etas(i,j-1)*((vy(i,j)-vy(i,j-1))*xkf+(vx(i+1,j-1)-vx(i,j-1))*xykf));
%                resy(i,j)=RY(i,j)-(xkf1*(vy(i,j-1)+vy(i,j+1))+ykf1*(vy(i-1,j)+vy(i+1,j))+xykf1*vy(i,j)-(pr(i,j-1)-pr(i-1,j-1))/ystp);

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


