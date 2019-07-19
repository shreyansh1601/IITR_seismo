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
% Layers thickness (m), density (kg/m3), viscosity (Pa s)
% Upper Layer
MH(1)=1500;
MRHO(1)=3000.0;
META(1)=1e+13;
% Lower Layer
MH(2)=1500;
MRHO(2)=2900.0;
META(2)=1e+19;

% Perturbation wavelength, amplitude, m
lamb=4000;
dApertr=MH(1)/150

% Model size, m
xsize=2*lamb;
ysize=MH(1)+MH(2);

% Velocity Boundary condition specified by bleft,bright,btop,bbot 
% (1=free slip -1=no slip) are implemented from ghost nodes 
% directly into Stokes and continuity equations
bleft=1;
bright=1;
btop=-1;
bbottom=-1;

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
        % Density, viscosity structure layering
        if(MY(ym,xm)>=MH(1))
            MI(ym,xm)=2;
        end
        % Cosinusoidal vertical perturbation
        MY(ym,xm)=MY(ym,xm)-cos(2*pi*(MX(ym,xm)-xsize/2)/lamb)*dApertr;
    end
end




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
rho1=rho1./wtnodes;
typ1=typ1./wtetas;
etas1=etas1./wtetas;
etan1=etan1./wtetan;

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
%     title(['x-Stokes residual, Gauss-Seidel cycle = ',num2str(niter)])
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
    surf(rho1);
    shading interp;
    light;
    lighting phong;
    axis tight;
%     colorbar;


% Define vertical velocity for the top of the growing wave
% Position of the top of the wave in vy grid
xwave=(xsize+xstp)/2;
ywave=MH(1)-dApertr;
% Define indexes for upper left node in the Vy cell where the top of the wave is
% !!! SUBTRACT 0.5 since int16(0.5)=1
xn=double(int16(xwave./xstp-0.5))+1;
yn=double(int16(ywave./ystp-0.5))+1;
% Define normalized distances from the top of the growing wave to the upper left Vy node;
dx=xwave./xstp-xn+1;
dy=ywave./ystp-yn+1;
% Interpolate Vy velocity from 4 nodes
wvy=(1.0-dx)*(1.0-dy)*vy1(yn,xn)+dx*(1.0-dy)*vy1(yn,xn+1)+(1.0-dx)*dy*vy1(yn+1,xn)+dx*dy*vy1(yn+1,xn+1)
Q=(MRHO(1)-MRHO(2))*MH(2)*g/2/META(2)
% Compute growth factor K
K=abs(wvy)/dApertr/Q


