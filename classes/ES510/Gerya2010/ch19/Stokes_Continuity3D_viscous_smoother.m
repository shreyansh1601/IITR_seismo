% Function Stokes_Continuity3D_viscous_smoother()
% This function makes specified number (iternum) of Gauss Seidel iteraions 
% using relaxation koefficients (krelaxs,krelaxc)
% for Stokes and Continuity equations defined on 2D staggered grid
% with specified resolution (xnum,ynum,znum) and gridsteps (xstp,ystp,zstp)
% given distribution of right parts for all equations (RX,RY,RZ,RC) on the grid 
% and given variable viscosity (etaxy,etaxz,etayz,etan) on the grid 
% pressure is normalized relative to given value (prnorm) in the first cell 
%
% Function return new approximation for velocity and pressure (vx,vy,vz,pr)
% and distribution of residuals (resx,resy,resz,resc)
function[vx,resx,vy,resy,vz,resz,pr,resc]=Stokes_Continuity3D_viscous_smoother(prnorm,etaxy,etaxz,etayz,etan,iternum,krelaxs,krelaxc,xnum,ynum,znum,xstp,ystp,zstp,RX,vx,RY,vy,RZ,vz,RC,pr)
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
ykf=1/ystp^2;
zkf=1/zstp^2;
xkf2=2/xstp^2;
ykf2=2/ystp^2;
zkf2=2/zstp^2;
xykf=1/xstp/ystp;
xzkf=1/xstp/zstp;
yzkf=1/ystp/zstp;

% Gauss-Seidel vx, xy, xz, P iteration circle
for niter=1:1:iternum;
    % Solving of Stokes and continuity equations on nodes
    for i=1:1:ynum+1;
        for j=1:1:xnum+1;
            for k=1:1:znum+1;

                % x-Stokes equation dSIGMAxx/dx+dSIGMAxy/dy+SIGMAxz/dz-dP/dx=RX
                if (j<xnum+1)
                    % vx-Boundrary conditions 
                    if (i==1 || i==ynum+1 || j==1 || j==xnum || k==1 || k==znum+1 );
                        % Boundary condition for vx
                        vx(i,j,k)=0;
                        % Upper Boundary
                        if(i==1)
                            % Free slip
                            vx(i,j,k)=vx(i+1,j,k);
                            % No slip
                            %vx(i,j,k)=-vx(i+1,j,k);
                        end
                        % Lower boundary
                        if(i==ynum+1)
                            % Free slip
                            vx(i,j,k)=vx(i-1,j,k);
                            % No slip
                            %vx(i,j,k)=-vx(i-1,j,k);
                        end
                        % Back Boundary
                        if(k==1)
                            % Free slip
                            vx(i,j,k)=vx(i,j,k+1);
                            % No slip
                            %vx(i,j,k)=-vx(i,j,k+1);
                        end
                        % Front boundary
                        if(k==znum+1)
                            % Free slip
                            vx(i,j,k)=vx(i,j,k-1);
                            % No slip
                            %vx(i,j,k)=-vx(i,j,k-1);
                        end
                    % Solving x-Stokes equation
                    else
                        % Computing current residual
                        resxcur=RX(i,j,k)+(pr(i-1,j,k-1)-pr(i-1,j-1,k-1))/xstp;
                        % dSIGMAxx/dx
                        resxcur=resxcur-(xkf2*(etan(i-1,j,k-1)*(vx(i,j+1,k)-vx(i,j,k))-etan(i-1,j-1,k-1)*(vx(i,j,k)-vx(i,j-1,k))));
                        % dSIGMAxy/dy
                        resxcur=resxcur-(etaxy(i,j,k-1)*(ykf*(vx(i+1,j,k)-vx(i,j,k))+xykf*(vy(i,j+1,k)-vy(i,j,k)))-etaxy(i-1,j,k-1)*(ykf*(vx(i,j,k)-vx(i-1,j,k))+xykf*(vy(i-1,j+1,k)-vy(i-1,j,k))));
                        % dSIGMAxz/dz
                        resxcur=resxcur-(etaxz(i-1,j,k)*(zkf*(vx(i,j,k+1)-vx(i,j,k))+xzkf*(vz(i,j+1,k)-vz(i,j,k)))-etaxz(i-1,j,k-1)*(zkf*(vx(i,j,k)-vx(i,j,k-1))+xzkf*(vz(i,j+1,k-1)-vz(i,j,k-1))));
                        % Koefficient at vx(i,j,k)
                        kfxcur=-xkf2*(etan(i-1,j,k-1)+etan(i-1,j-1,k-1))-ykf*(etaxy(i,j,k-1)+etaxy(i-1,j,k-1))-zkf*(etaxz(i-1,j,k)+etaxz(i-1,j,k-1));
                        % Updating solution
                        vx(i,j,k)=vx(i,j,k)+resxcur/kfxcur*krelaxs;

                    end
                end
            
                % y-Stokes equation dSIGMAyx/dx+dSIGMAyy/dy+SIGMAyz/dz-dP/dy=RY
                if (i<ynum+1)
                    % vy-Boundrary conditions 
                    if (i==1 || i==ynum || j==1 || j==xnum+1 || k==1 || k==znum+1);
                    % Boundary condition for vy
                        vy(i,j,k)=0;
                        % Left boundary
                        if(j==1)
                            % Free slip
                            vy(i,j,k)=vy(i,j+1,k);
                            % No slip
                            %vy(i,j,k)=-vy(i,j+1,k);
                        end
                        % Right boundary
                        if(j==xnum+1)
                            % Free slip
                            vy(i,j,k)=vy(i,j-1,k);
                            % No slip
                            %vy(i,j,k)=-vy(i,j-1,k);
                        end
                        % Back Boundary
                        if(k==1)
                            % Free slip
                            vy(i,j,k)=vy(i,j,k+1);
                            % No slip
                            %vy(i,j,k)=-vy(i,j,k+1);
                        end
                        % Front boundary
                        if(k==znum+1)
                            % Free slip
                            vy(i,j,k)=vy(i,j,k-1);
                            % No slip
                            %vy(i,j,k)=-vy(i,j,k-1);
                        end
                    % Solving y-Stokes equation
                    else
                        % Computing current residual
                        resycur=RY(i,j,k)+(pr(i,j-1,k-1)-pr(i-1,j-1,k-1))/ystp;
                        % dSIGMAyy/dy
                        resycur=resycur-(ykf2*(etan(i,j-1,k-1)*(vy(i+1,j,k)-vy(i,j,k))-etan(i-1,j-1,k-1)*(vy(i,j,k)-vy(i-1,j,k))));
                        % dSIGMAyx/dx
                        resycur=resycur-(etaxy(i,j,k-1)*(xkf*(vy(i,j+1,k)-vy(i,j,k))+xykf*(vx(i+1,j,k)-vx(i,j,k)))-etaxy(i,j-1,k-1)*(xkf*(vy(i,j,k)-vy(i,j-1,k))+xykf*(vx(i+1,j-1,k)-vx(i,j-1,k))));
                        % dSIGMAyz/dz
                        resycur=resycur-(etayz(i,j-1,k)*(zkf*(vy(i,j,k+1)-vy(i,j,k))+yzkf*(vz(i+1,j,k)-vz(i,j,k)))-etayz(i,j-1,k-1)*(zkf*(vy(i,j,k)-vy(i,j,k-1))+yzkf*(vz(i+1,j,k-1)-vz(i,j,k-1))));
                        % Koefficient at vy(i,j,k)
                        kfycur=-ykf2*(etan(i,j-1,k-1)+etan(i-1,j-1,k-1))-xkf*(etaxy(i,j,k-1)+etaxy(i,j-1,k-1))-zkf*(etayz(i,j-1,k)+etayz(i,j-1,k-1));
                        % Updating solution
                        vy(i,j,k)=vy(i,j,k)+resycur/kfycur*krelaxs;
                        
                    end
                end
            
                % z-Stokes equation dSIGMAzx/dx+dSIGMAzy/dy+SIGMAzz/dz-dP/dz=RZ
                if (k<znum+1)
                    % vz-Boundrary conditions 
                    if (i==1 || i==ynum+1 || j==1 || j==xnum+1 || k==1 || k==znum);
                    % Boundary condition for vz
                        vz(i,j,k)=0;
                        % Upper Boundary
                        if(i==1)
                            % Free slip
                            vz(i,j,k)=vz(i+1,j,k);
                            % No slip
                            %vz(i,j,k)=-vz(i+1,j,k);
                        end
                        % Lower boundary
                        if(i==ynum+1)
                            % Free slip
                            vz(i,j,k)=vz(i-1,j,k);
                            % No slip
                            %vz(i,j,k)=-vz(i-1,j,k);
                        end
                        % Left boundary
                        if(j==1)
                            % Free slip
                            vz(i,j,k)=vz(i,j+1,k);
                            % No slip
                            %vz(i,j,k)=-vz(i,j+1,k);
                        end
                        % Right boundary
                        if(j==xnum+1)
                            % Free slip
                            vz(i,j,k)=vz(i,j-1,k);
                            % No slip
                            %vz(i,j,k)=-vz(i,j-1,k);
                        end
                    % Solving z-Stokes equation
                    else
                        % Computing current residual
                        reszcur=RZ(i,j,k)+(pr(i-1,j-1,k)-pr(i-1,j-1,k-1))/zstp;
                        % dSIGMAzz/dz
                        reszcur=reszcur-(zkf2*(etan(i-1,j-1,k)*(vz(i,j,k+1)-vz(i,j,k))-etan(i-1,j-1,k-1)*(vz(i,j,k)-vz(i,j,k-1))));
                        % dSIGMAzx/dx
                        reszcur=reszcur-(etaxz(i-1,j,k)*(xkf*(vz(i,j+1,k)-vz(i,j,k))+xzkf*(vx(i,j,k+1)-vx(i,j,k)))-etaxz(i-1,j-1,k)*(xkf*(vz(i,j,k)-vz(i,j-1,k))+xzkf*(vx(i,j-1,k+1)-vx(i,j-1,k))));
                        % dSIGMAzy/dy
                        reszcur=reszcur-(etayz(i,j-1,k)*(ykf*(vz(i+1,j,k)-vz(i,j,k))+yzkf*(vy(i,j,k+1)-vy(i,j,k)))-etayz(i-1,j-1,k)*(ykf*(vz(i,j,k)-vz(i-1,j,k))+yzkf*(vy(i-1,j,k+1)-vy(i-1,j,k))));
                        % Koefficient at vy(i,j,k)
                        kfzcur=-zkf2*(etan(i-1,j-1,k)+etan(i-1,j-1,k-1))-xkf*(etaxz(i-1,j,k)+etaxz(i-1,j-1,k))-ykf*(etayz(i,j-1,k)+etayz(i-1,j-1,k));
                        % Updating solution
                        vz(i,j,k)=vz(i,j,k)+reszcur/kfzcur*krelaxs;
                    end
                end
                
                % Continuity equation dvx/dx+dvy/dy+dvz/dz=RC
                % is solved via pressure updates
                % dpr=-etas*div(v)
                % pr-Boundrary conditions 
                if (i<ynum && j<xnum && k<znum);
                    % Solving Continuity equation by adjusting pressure
                    % Computing current residual
                    resc(i,j,k)=RC(i,j,k)-((vx(i+1,j+1,k+1)-vx(i+1,j,k+1))/xstp+(vy(i+1,j+1,k+1)-vy(i,j+1,k+1))/ystp+(vz(i+1,j+1,k+1)-vz(i+1,j+1,k))/zstp);
                    % Updating pressure solution
                    pr(i,j,k)=pr(i,j,k)+resc(i,j,k)*etan(i,j,k)*krelaxc;
                end
                
            end
        end            
    end
    % End of Solving Stokes and continuity equations on nodes
end

% Correct pressure with the given value in the first cell
dp=prnorm-pr(1,1,1);
pr=pr+dp;

% Computing final state of residuals
for i=1:1:ynum+1;
    for j=1:1:xnum+1;
        for k=1:1:znum+1;

            % x-Stokes equation d2vx/dx2+d2vx/dy2+d2vx/dz2-dP/dx=RX
            if (j<xnum+1)
                % vx-Boundrary conditions 
                if (i==1 || i==ynum+1 || j==1 || j==xnum || k==1 || k==znum+1 );
                    resx(i,j,k)=0;
                % x-Stokes equation
                else
                    % Computing current residual
                    resxcur=RX(i,j,k)+(pr(i-1,j,k-1)-pr(i-1,j-1,k-1))/xstp;
                    % dSIGMAxx/dx
                    resxcur=resxcur-(xkf2*(etan(i-1,j,k-1)*(vx(i,j+1,k)-vx(i,j,k))-etan(i-1,j-1,k-1)*(vx(i,j,k)-vx(i,j-1,k))));
                    % dSIGMAxy/dy
                    resxcur=resxcur-(etaxy(i,j,k-1)*(ykf*(vx(i+1,j,k)-vx(i,j,k))+xykf*(vy(i,j+1,k)-vy(i,j,k)))-etaxy(i-1,j,k-1)*(ykf*(vx(i,j,k)-vx(i-1,j,k))+xykf*(vy(i-1,j+1,k)-vy(i-1,j,k))));
                    % dSIGMAxz/dz
                    resxcur=resxcur-(etaxz(i-1,j,k)*(zkf*(vx(i,j,k+1)-vx(i,j,k))+xzkf*(vz(i,j+1,k)-vz(i,j,k)))-etaxz(i-1,j,k-1)*(zkf*(vx(i,j,k)-vx(i,j,k-1))+xzkf*(vz(i,j+1,k-1)-vz(i,j,k-1))));
                    resx(i,j,k)=resxcur;
                end
            end
            
            % y-Stokes equation d2vy/dx2+d2vy/dy2+d2vy/dz2-dP/dy=RY
            if (i<ynum+1)
                % vy-Boundrary conditions 
                if (i==1 || i==ynum || j==1 || j==xnum+1 || k==1 || k==znum+1);
                resy(i,j,k)=0;
                %y-Stokes equation
                else
                    % Computing current residual
                    resycur=RY(i,j,k)+(pr(i,j-1,k-1)-pr(i-1,j-1,k-1))/ystp;
                    % dSIGMAyy/dy
                    resycur=resycur-(ykf2*(etan(i,j-1,k-1)*(vy(i+1,j,k)-vy(i,j,k))-etan(i-1,j-1,k-1)*(vy(i,j,k)-vy(i-1,j,k))));
                    % dSIGMAyx/dx
                    resycur=resycur-(etaxy(i,j,k-1)*(xkf*(vy(i,j+1,k)-vy(i,j,k))+xykf*(vx(i+1,j,k)-vx(i,j,k)))-etaxy(i,j-1,k-1)*(xkf*(vy(i,j,k)-vy(i,j-1,k))+xykf*(vx(i+1,j-1,k)-vx(i,j-1,k))));
                    % dSIGMAyz/dz
                    resycur=resycur-(etayz(i,j-1,k)*(zkf*(vy(i,j,k+1)-vy(i,j,k))+yzkf*(vz(i+1,j,k)-vz(i,j,k)))-etayz(i,j-1,k-1)*(zkf*(vy(i,j,k)-vy(i,j,k-1))+yzkf*(vz(i+1,j,k-1)-vz(i,j,k-1))));
                    resy(i,j,k)=resycur;
                end
            end
            
            % z-Stokes equation d2vz/dx2+d2vz/dy2+d2vz/dz2-dP/dz=RZ
            if (k<znum+1)
                % vz-Boundrary conditions 
                if (i==1 || i==ynum+1 || j==1 || j==xnum+1 || k==1 || k==znum);
                    resz(i,j,k)=0;
                %y-Stokes equation
                else
                    % Computing current residual
                    reszcur=RZ(i,j,k)+(pr(i-1,j-1,k)-pr(i-1,j-1,k-1))/zstp;
                    % dSIGMAzz/dz
                    reszcur=reszcur-(zkf2*(etan(i-1,j-1,k)*(vz(i,j,k+1)-vz(i,j,k))-etan(i-1,j-1,k-1)*(vz(i,j,k)-vz(i,j,k-1))));
                    % dSIGMAzx/dx
                    reszcur=reszcur-(etaxz(i-1,j,k)*(xkf*(vz(i,j+1,k)-vz(i,j,k))+xzkf*(vx(i,j,k+1)-vx(i,j,k)))-etaxz(i-1,j-1,k)*(xkf*(vz(i,j,k)-vz(i,j-1,k))+xzkf*(vx(i,j-1,k+1)-vx(i,j-1,k))));
                    % dSIGMAzy/dy
                    reszcur=reszcur-(etayz(i,j-1,k)*(ykf*(vz(i+1,j,k)-vz(i,j,k))+yzkf*(vy(i,j,k+1)-vy(i,j,k)))-etayz(i-1,j-1,k)*(ykf*(vz(i,j,k)-vz(i-1,j,k))+yzkf*(vy(i-1,j,k+1)-vy(i-1,j,k))));
                    resz(i,j,k)=reszcur;
                end
            end

            % Continuity equation dvx/dx+dvy/dy+dvz/dz=RC
            % is solved via pressure updates
            % dpr=-etas*div(v)
            if (i<ynum && j<xnum && k<znum);
                % Computing current residual
                resc(i,j,k)=RC(i,j,k)-((vx(i+1,j+1,k+1)-vx(i+1,j,k+1))/xstp+(vy(i+1,j+1,k+1)-vy(i,j+1,k+1))/ystp+(vz(i+1,j+1,k+1)-vz(i+1,j+1,k))/zstp);
            end
        
        end
    end            
end


