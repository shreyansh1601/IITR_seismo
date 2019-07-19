% Solution of Stokes and continuity equations 
% with variable viscosity in 2D with direct solver
% by using external function Stokes_Continuity_solver_ghost()
% Setup corresponds to falling block test
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

% Clearing all variables and arrays
clear all;
% Clearing figures
clf;


% Acceleration of Gravity, m/s^2
g=9.81;
% Pressure in the upermost, leftmost (first) cell
prfirst=0;
% Rock density (kg/m3), viscosity (Pa s)
% Medium
MRHO(1)=3200.0;
META(1)=1e+21;
% Block
MRHO(2)=3300.0;
META(2)=1e+27;

% Maximal timestep, s
timemax=1e+8*(365.25*24*3600);
% Maximal marker displacement step, number of gridsteps
markmax=0.5;
% Amount of timesteps
stepmax=100;

% Model size, m
xsize=500000;
ysize=500000;

% Velocity Boundary condition specified by bleft,bright,btop,bbot 
% (1=free slip -1=no slip) are implemented from ghost nodes 
% directly into Stokes and continuity equations
bleft=1;
bright=1;
btop=1;
bbottom=1;

% Defining resolution
xnum=51;
ynum=51;

% Defining gridsteps
xstp=xsize./(xnum-1);
ystp=ysize./(ynum-1);




% Defining number of markers and steps between them in the horizontal and vertical direction
xmx=5; %number of markers per cell in horizontal direction
ymy=5; %number of markers per cell in vertical direction
mxnum=(xnum-1)*xmx; %total number of markers in horizontal direction
mynum=(ynum-1)*ymy; %total number of markers in vertical direction
mxstep=xsize/mxnum; %step between markers in horizontal direction   
mystep=ysize/mynum; %step between markers in vertical direction

% Creating markers arrays
MX=zeros(mynum,mxnum); % X coordinate
MY=zeros(mynum,mxnum); % Y coordinate
MI=zeros(mynum,mxnum); % Type

% Defining intial position of markers
% Defining lithological structure of the model
for xm = 1:1:mxnum
    for ym = 1:1:mynum  
        MX(ym,xm)=xm*mxstep-mxstep/2;
        MY(ym,xm)=ym*mystep-mystep/2;
        MI(ym,xm)=1;
        % Density, viscosity structure definition for block
        % Relative distances for the marker inside the grid
        dx=MX(ym,xm)/xsize;
        dy=MY(ym,xm)/ysize;
        if(dx>=0.4 && dx<=0.6 && dy>=0.1 && dy<=0.3)
            MI(ym,xm)=2;
        end
    end
end


% Rock type, density and viscosity arrays
typ1 = zeros(ynum,xnum);
etas1 = zeros(ynum,xnum);
etan1 = zeros(ynum-1,xnum-1);
rho1 = zeros(ynum,xnum);

% Initial time, s
timesum=0;

% Main Time cycle
% Backup rock type, density and viscosity arrays
for ntimestep=1:1:stepmax
    typ0 = typ1;
    etas0 = etas1;
    etan0 = etan1;
    rho0 = rho1;
    % Clear rock type, density and viscosity arrays
    typ1 = zeros(ynum,xnum);
    etas1 = zeros(ynum,xnum);
    etan1 = zeros(ynum-1,xnum-1);
    rho1 = zeros(ynum,xnum);
    % Clear wights for basic nodes
    wtnodes=zeros(ynum,xnum);
    % Clear wights for etas
    wtetas=zeros(ynum,xnum);
    % Clear wights for etan
    wtetan=zeros(ynum-1,xnum-1);

    % Interpolating parameters from markers to nodes
    for xm = 1:1:mxnum
        for ym = 1:1:mynum  

            %  xn    rho(xn,yn)--------------------rho(xn+1,yn)
            %           ?           ^                  ?
            %           ?           ?                  ?
            %           ?          dy                  ?
            %           ?           ?                  ?
            %           ?           v                  ?
            %           ?<----dx--->o Mrho(xm,ym)       ?
            %           ?                              ?
            %           ?                              ?
            %  xn+1  rho(xn,yn+1)-------------------rho(xn+1,yn+1)
            %
            % Define indexes for upper left node in the cell where the marker is
            % !!! SUBTRACT 0.5 since int16(0.5)=1
            xn=double(int16(MX(ym,xm)./xstp-0.5))+1;
            yn=double(int16(MY(ym,xm)./ystp-0.5))+1;
            if (xn<1)
                xn=1;
            end
            if (xn>xnum-1)
                xn=xnum-1;
            end
            if (yn<1)
                yn=1;
            end
            if (yn>ynum-1)
                yn=ynum-1;
            end

            % Define normalized distances from marker to the upper left node;
            dx=MX(ym,xm)./xstp-xn+1;
            dy=MY(ym,xm)./ystp-yn+1;

            % Add density to 4 surrounding nodes
            rho1(yn,xn)=rho1(yn,xn)+(1.0-dx).*(1.0-dy).*MRHO(MI(ym,xm));
            wtnodes(yn,xn)=wtnodes(yn,xn)+(1.0-dx).*(1.0-dy);
            rho1(yn+1,xn)=rho1(yn+1,xn)+(1.0-dx).*dy.*MRHO(MI(ym,xm));
            wtnodes(yn+1,xn)=wtnodes(yn+1,xn)+(1.0-dx).*dy;
            rho1(yn,xn+1)=rho1(yn,xn+1)+dx.*(1.0-dy).*MRHO(MI(ym,xm));
            wtnodes(yn,xn+1)=wtnodes(yn,xn+1)+dx.*(1.0-dy);
            rho1(yn+1,xn+1)=rho1(yn+1,xn+1)+dx.*dy.*MRHO(MI(ym,xm));
            wtnodes(yn+1,xn+1)=wtnodes(yn+1,xn+1)+dx.*dy;

            % Add shear viscosity etas() and rock type typ() to 4 surrounding nodes
            % only using markers located at <=0.5 gridstep distances from nodes
            if(dx<=0.5 && dy<=0.5)
                etas1(yn,xn)=etas1(yn,xn)+(1.0-dx).*(1.0-dy).*META(MI(ym,xm));
                typ1(yn,xn)=typ1(yn,xn)+(1.0-dx).*(1.0-dy).*MI(ym,xm);
                wtetas(yn,xn)=wtetas(yn,xn)+(1.0-dx).*(1.0-dy);
            end
            if(dx<=0.5 && dy>=0.5)
                etas1(yn+1,xn)=etas1(yn+1,xn)+(1.0-dx).*dy.*META(MI(ym,xm));
                typ1(yn+1,xn)=typ1(yn+1,xn)+(1.0-dx).*dy.*MI(ym,xm);
                wtetas(yn+1,xn)=wtetas(yn+1,xn)+(1.0-dx).*dy;
            end
            if(dx>=0.5 && dy<=0.5)
                etas1(yn,xn+1)=etas1(yn,xn+1)+dx.*(1.0-dy).*META(MI(ym,xm));
                typ1(yn,xn+1)=typ1(yn,xn+1)+dx.*(1.0-dy).*MI(ym,xm);
                wtetas(yn,xn+1)=wtetas(yn,xn+1)+dx.*(1.0-dy);
            end
            if(dx>=0.5 && dy>=0.5)
                etas1(yn+1,xn+1)=etas1(yn+1,xn+1)+dx.*dy.*META(MI(ym,xm));
                typ1(yn+1,xn+1)=typ1(yn+1,xn+1)+dx.*dy.*MI(ym,xm);
                wtetas(yn+1,xn+1)=wtetas(yn+1,xn+1)+dx.*dy;
            end

            % Add normal viscosity etan() to the center of current cell
            etan1(yn,xn)=etan1(yn,xn)+(1.0-abs(0.5-dx)).*(1.0-abs(0.5-dy)).*META(MI(ym,xm));
            wtetan(yn,xn)=wtetan(yn,xn)+(1.0-abs(0.5-dx)).*(1.0-abs(0.5-dy));

        end
    end

    % Computing  Viscosity, density, rock type for nodal points
    for i=1:1:ynum;
        for j=1:1:xnum;
            % Density
            if (wtnodes(i,j)~=0)
                % Compute new value interpolated from markers
                rho1(i,j)=rho1(i,j)./wtnodes(i,j);
            else
                % If no new value is interpolated from markers old value is used
                rho1(i,j)=rho0(i,j);
            end
            % Shear viscosity and type
            if (wtetas(i,j)~=0)
                % Compute new value interpolated from markers
                etas1(i,j)=etas1(i,j)./wtetas(i,j);
                typ1(i,j)=typ1(i,j)./wtetas(i,j);
            else
                % If no new value is interpolated from markers old value is used
                etas1(i,j)=etas0(i,j);
                typ1(i,j)=typ0(i,j);
            end
            % Normal viscosity
            if (i<ynum && j<xnum)
                if (wtetan(i,j)~=0)
                    % Compute new value interpolated from markers
                    etan1(i,j)=etan1(i,j)./wtetan(i,j);
                else
                    % If no new value is interpolated from markers old value is used
                    etan1(i,j)=etan0(i,j);
                end
            end
        end
    end


    % Computing right part of Stokes (RX, RY) and Continuity (RC) equation
    % vx, vy, P
    vx1=zeros(ynum+1,xnum);
    vy1=zeros(ynum,xnum+1);
    pr1=zeros(ynum-1,xnum-1);
    % Right parts of equations
    RX1=zeros(ynum+1,xnum);
    RY1=zeros(ynum,xnum+1);
    RC1=zeros(ynum-1,xnum-1);
    % Grid points cycle
    for i=1:1:ynum;
        for j=1:1:xnum;
            % Right part of x-Stokes Equation
            if(j>1 && i>1 && j<xnum)
                RX1(i,j)=0;
            end
            % Right part of y-Stokes Equation
            if(j>1 && i>1 && i<ynum)
                RY1(i,j)=-g*(rho1(i,j)+rho1(i,j-1))/2;
            end
        end
    end





    % Solving of Stokes and Continuity equations on nodes
    % and computing residuals
    % by calling function Stokes_Continuity_solver_ghost()

    [vx1,resx1,vy1,resy1,pr1,resc1]=Stokes_Continuity_solver_ghost(prfirst,etas1,etan1,xnum,ynum,xstp,ystp,RX1,RY1,RC1,bleft,bright,btop,bbottom);

    % Defining scale for Stokes residuals from y-Stokes equation
    % dSIGMAij/dj-dP/di=-RHO*gi=0  => Stokes scale=abs(RHO*gi)
    stokesscale= MRHO(1)*g;
    % Defining scale for Continuity residuals from y-Stokes equation
    % dvx/dx+dvy/dy=0 can be transformed to 2ETA(dvx/dx+dvy/dy)/dx=0 
    % which is similar to dSIGMAij/dj and has scale given above
    % therefore continuity scale = scale=abs(RHO*gi/ETA*dx)
    continscale= MRHO(1)*g/META(1)*xstp;



    % Defining timestep
    timestep=timemax;
    % Check maximal velocity
    vxmax=max(abs(max(max(vx1))),abs(min(min(vx1))))
    vymax=max(abs(max(max(vy1))),abs(min(min(vy1))))
    % Check marker displacement step
    if (vxmax>0)
        if (timestep>markmax*xstp/vxmax)
            timestep=markmax*xstp/vxmax
        end
    end
    if (vymax>0)
        if (timestep>markmax*ystp/vymax)
            timestep=markmax*ystp/vymax
        end
    end
    timestep=timestep

    % Moving Markers by velocity field
    for xm = 1:1:mxnum
        for ym = 1:1:mynum  

            %  xn    V(xn,yn)--------------------V(xn+1,yn)
            %           ?           ^                  ?
            %           ?           ?                  ?
            %           ?          dy                  ?
            %           ?           ?                  ?
            %           ?           v                  ?
            %           ?<----dx--->o Mrho(xm,ym)       ?
            %           ?                              ?
            %           ?                              ?
            %  xn+1  V(xn,yn+1)-------------------V(xn+1,yn+1)
            %
            % Define indexes for upper left node in the VX-cell where the marker is
            % VX-cells are displaced upward for 1/2 of vertical gridstep
            % !!! SUBTRACT 0.5 since int16(0.5)=1
            xn=double(int16(MX(ym,xm)./xstp-0.5))+1;
            yn=double(int16((MY(ym,xm)+ystp/2.0)./ystp-0.5))+1;
            % Check vertical index for upper left VX-node 
            % It must be between 1 and ynum (see picture for staggered grid)
            if (xn<1)
                xn=1;
            end
            if (xn>xnum-1)
                xn=xnum-1;
            end
            if (yn<1)
                yn=1;
            end
            if (yn>ynum)
                yn=ynum;
            end
            % Define and check normalized distances from marker to the upper left VX-node;
            dx=MX(ym,xm)./xstp-xn+1;
            dy=(MY(ym,xm)+ystp/2.0)./ystp-yn+1;

            % Calculate Marker velocity from four surrounding nodes
            vxm=0;
            vxm=vxm+(1.0-dx).*(1.0-dy).*vx1(yn,xn);
            vxm=vxm+(1.0-dx).*dy.*vx1(yn+1,xn);
            vxm=vxm+dx.*(1.0-dy).*vx1(yn,xn+1);
            vxm=vxm+dx.*dy.*vx1(yn+1,xn+1);

            % Define indexes for upper left node in the VY-cell where the marker is
            % VY-cells are displaced leftward for 1/2 of horizontal gridstep
            % !!! SUBTRACT 0.5 since int16(0.5)=1
            xn=double(int16((MX(ym,xm)+xstp/2.0)./xstp-0.5))+1;
            yn=double(int16(MY(ym,xm)./ystp-0.5))+1;
            % Check horizontal index for upper left VY-node 
            % It must be between 1 and xnum (see picture for staggered grid)
            if (xn<1)
                xn=1;
            end
            if (xn>xnum)
                xn=xnum;
            end
            if (yn<1)
                yn=1;
            end
            if (yn>ynum-1)
                yn=ynum-1;
            end
            % Define and check normalized distances from marker to the upper left VX-node;
            dx=(MX(ym,xm)+xstp/2.0)./xstp-xn+1;
            dy=MY(ym,xm)./ystp-yn+1;
            % Calculate Marker velocity from four surrounding nodes
            vym=0;
            vym=vym+(1.0-dx).*(1.0-dy).*vy1(yn,xn);
            vym=vym+(1.0-dx).*dy.*vy1(yn+1,xn);
            vym=vym+dx.*(1.0-dy).*vy1(yn,xn+1);
            vym=vym+dx.*dy.*vy1(yn+1,xn+1);

            % Displacing Marker according to its velocity
            MX(ym,xm)=MX(ym,xm)+timestep*vxm;
            MY(ym,xm)=MY(ym,xm)+timestep*vym;
        end
    end
    % Plotting Residuals for x-Stokes as surface
    figure(1);
    subplot(2,3,1)
    resx0=resx1/stokesscale;
    surf(resx0);
    shading interp;
    light;
    lighting phong;
    axis tight;
    zlabel('residual x-Stokes')
%     colorbar;
    
    % Plotting Residuals for y-Stokes as surface
    subplot(2,3,2)
    resy0=resy1/stokesscale;
    surf(resy0);
    shading interp;
    light;
    lighting phong;
    axis tight;
    zlabel('residual Y-stokes')
%     colorbar;
    
    % Plotting Residuals for Continuity as surface
    subplot(2,3,3)
    resc0=resc1/continscale;
    surf(resc0);
    shading interp;
    light;
    lighting phong;
    axis tight;
    zlabel('residual continuity')
%     colorbar;
    
    % Plotting vx
    subplot(2,3,4)
    surf(vx1);
    shading interp;
    light;
    lighting phong;
    axis tight;
%     colorbar;

    % Plotting vy
    subplot(2,3,5)
    surf(vy1);
    shading interp;
    light;
    lighting phong;
    axis tight;
%     colorbar;

    % Plotting P
    subplot(2,3,6)
    pcolor(typ1);
    shading interp;
    colorbar;
    axis ij;
    title(['Rock types, Time = ',num2str(timesum*1e-6/(365.25*24*3600)),' Myr']);

    % Advance in time
    timesum=timesum+timestep;

    % Define vertical velocity for the initial center of the block
    % Initial position of block center in vy grid
    xblock=(xsize+xstp)/2;
    yblock=0.2*ysize;
    % Define indexes for upper left node in the Vy cell where the top of the wave is
    % !!! SUBTRACT 0.5 since int16(0.5)=1
    xn=double(int16(xblock./xstp-0.5))+1;
    yn=double(int16(yblock./ystp-0.5))+1;
    % Define normalized distances from the top of the growing wave to the upper left Vy node;
    dx=xblock./xstp-xn+1;
    dy=yblock./ystp-yn+1;
    % Interpolate Vy velocity from 4 nodes
    vyblock=(1.0-dx)*(1.0-dy)*vy1(yn,xn)+dx*(1.0-dy)*vy1(yn,xn+1)+(1.0-dx)*dy*vy1(yn+1,xn)+dx*dy*vy1(yn+1,xn+1)

end