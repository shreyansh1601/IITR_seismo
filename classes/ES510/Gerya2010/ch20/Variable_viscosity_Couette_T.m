% Solution of Stokes and continuity equations 
% with variable viscosity in 2D with direct solver
% by using external function Stokes_Continuity_solver_channel()
% solving temperature equation
% Setup corresponds to constant viscosity channel test
% with non steady temperature distribution
% 
% Staggered Grid 
% 
%     vx       vx       vx    
%
% vy  T---vy---T---vy---T   vy
%     |        |        |
%     vx   P   vx   P   vx    
%     |        |        |
% vy  T---vy---T---vy---T   vy
%     |        |        |
%     vx   P   vx   P   vx    
%     |        |        |
% vy  T---vy---T---vy---T   vy
%
%     vx       vx       vx    
% 
% Lines show basic grid
% Ghost nodes shown outside the basic grid
% are used for boundary conditions

% Clearing all variables and arrays
clear all;
% Clearing figures
clf;

% Temperature at the top, and bottom of the channel, K
ttop=1000;
tbottom=1000;

% Acceleration of Gravity, m/s^2
gx=0;
gy=0;

% Gas constant J/mol/K
RGAS=8.314;

% Pressure boundary conditions
% prfirst(1) = boundary condition mode:
% 0 - pressure in one cell definition
% 1 - pressure at the top and in the bottom of the channel
prfirst(1)=1;
% prfirst(2) = boundary condition value
prfirst(2)=0;

% Channel Medium
% Rock density (kg/m3)
MRHO(1)=3000.0;
% Rock flow law: 
% 0 = constant viscosity;
% 1 = isothermal power-law 2*EPSILONii=C1*SIGMAii^n
% 2 = temperature dependent ETA=A*exp(Ea/RT0(1-(T-T0)/T0))
MLAW(1)=2;
META(1,1)=1e+15; % A, Pa s
META(1,2)=1.5e+5; % Ea, J/mol
META(1,3)=ttop; % T0, K
% Rock heat capacity, J/K/kg
MCP(1)=1000;
% Rock thermal conductivity, W/m/K
MKT(1)=2;
% Maximal timestep, s
timemax=5e+6*(365.25*24*3600);
% Maximal marker displacement step, number of gridsteps
markmax=5;
% Moving Markers y(1)/n(0)
markmove=0;
% Maximal temperature change, allowed for one timestep, K
tempmax=10;
% Amount of timesteps
stepmax=100;

% Model size, m
xsize=30000;
ysize=12500;

% Defining resolution
xnum=51;
ynum=11;

% Velocity at the right wall
vyright=1e-9;
% Numerical Subgrid temperature diffusion coefficient
dsubgridt=0;

% Velocity Boundary condition specified by bleft,bright,btop,bbot 
% are implemented from ghost nodes 
% directly into Stokes and continuity equations
% Thermal boundary conditions
% Upper, Lower boundaries: symmetry
for j=1:1:xnum+1
    % Upper boundary:
    % vx(1,j)=btop(j,1)+vx(2,j)*btop(j,2)
    btop(j,1)=0;
    btop(j,2)=1;
    % vy(1,j)=btop(j,3)+vy(2,j)*btop(j,4)
    btop(j,3)=0;
    btop(j,4)=1;
    % Lower boundary
    % vx(ynum+1,j)=bbottom(j,1)+vx(ynum,j)*bbottom(j,2)
    bbottom(j,1)=0;
    bbottom(j,2)=1;
    % vy(ynum+,j)=bbottom(j,3)+vy(ynum-1,j)*bbottom(j,4)
    bbottom(j,3)=0;
    bbottom(j,4)=1;
end
% Left, Right boundaries: constant vertical velocity
for i=1:1:ynum+1
    % Left boundary: zero velocity
    % vx(i,1)=bleft(i,1)+vx(i,2)*bleft(i,2)
    bleft(i,1)=0;
    bleft(i,2)=0;
    % vy(i,1)=bleft(i,3)+vy(i,2)*bleft(i,42)
    bleft(i,3)=0;
    bleft(i,4)=-1;
    % Right boundary: constant vwrtical velocity
    % vx(i,xnum)=bright(i,1)+vx(i,xnum-1)*bbright(i,2)
    bright(i,1)=0;
    bright(i,2)=0;
    % vy(i,xnum+1)=bright(i,3)+vx(i,xnum)*bbright(i,4)
    bright(i,3)=2*vyright;
    bright(i,4)=-1;
end

% Defining gridsteps
xstp=xsize./(xnum-1);
ystp=ysize./(ynum-1);


% Vertical thermal gradient, K/m
dtdy=(tbottom-ttop)/ysize;

% Thermal boundary conditions
% Upper, Lower boundaries: constant thermal gradient
for j=1:1:xnum
    % Upper boundary
    btopt(j,1)=-ystp*dtdy;
    btopt(j,2)=1;
    % Lower boundary
    bbottomt(j,1)=ystp*dtdy;
    bbottomt(j,2)=1;
end
% Left, Right boundaries: constant temperature, symmetry 
for i=1:1:ynum
    % Left boundary: T= constant
    bleftt(i,1)=ttop+(i-1)*ystp*dtdy;
    bleftt(i,2)=0;
    % Right boundary: Symmetry
    brightt(i,1)=0;
    brightt(i,2)=1;
end


% Defining number of markers and steps between them in the horizontal and vertical direction
xmx=5; %number of markers per cell in horizontal direction
ymy=5; %number of markers per cell in vertical direction
mxnum=(xnum-1)*xmx; %total number of markers in horizontal direction
mynum=(ynum-1)*ymy; %total number of markers in vertical direction
mxstep=xsize/mxnum; %step between markers in horizontal direction   
mystep=ysize/mynum; %step between markers in vertical direction

% Defining intial position of markers
% Defining lithological structure of the model
for xm = 1:1:mxnum
    for ym = 1:1:mynum
        % Coordinates
        MX(ym,xm)=xm*mxstep-mxstep/2;
        MY(ym,xm)=ym*mystep-mystep/2;
        % Material Type
        MI(ym,xm)=1;
        % Second invariant of strain rate
        MEII(ym,xm)=0;
        % Initial temperature
        MTK(ym,xm)=ttop+MY(ym,xm)*dtdy;
    end
end


% Rock type, density, viscosity, temperature, thermal conductivity and RHO*Cp arrays
typ1 = zeros(ynum,xnum);
etas1 = zeros(ynum,xnum);
etan1 = zeros(ynum-1,xnum-1);
rho1 = zeros(ynum,xnum);
tk1 = zeros(ynum,xnum);
rhocp1 = zeros(ynum,xnum);
kt1 = zeros(ynum,xnum);

% Initial time, s
timesum=0;

% Main Time cycle
for ntimestep=1:1:stepmax


    % Backup transport properties arrays
    typ0 = typ1;
    etas0 = etas1;
    etan0 = etan1;
    rho0 = rho1;
    tk0=tk1;
    rhocp0=rhocp1;
    kt0=kt1;
    % Clear transport properties arrays
    typ1 = zeros(ynum,xnum);
    etas1 = zeros(ynum,xnum);
    etan1 = zeros(ynum-1,xnum-1);
    rho1 = zeros(ynum,xnum);
    tk1 = zeros(ynum,xnum);
    rhocp1 = zeros(ynum,xnum);
    kt1 = zeros(ynum,xnum);
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
            tk1(yn,xn)=tk1(yn,xn)+(1.0-dx).*(1.0-dy).*MTK(ym,xm);
            kt1(yn,xn)=kt1(yn,xn)+(1.0-dx).*(1.0-dy).*MKT(MI(ym,xm));
            rhocp1(yn,xn)=rhocp1(yn,xn)+(1.0-dx).*(1.0-dy).*MRHO(MI(ym,xm))*MCP(MI(ym,xm));
            wtnodes(yn,xn)=wtnodes(yn,xn)+(1.0-dx).*(1.0-dy);
            rho1(yn+1,xn)=rho1(yn+1,xn)+(1.0-dx).*dy.*MRHO(MI(ym,xm));
            tk1(yn+1,xn)=tk1(yn+1,xn)+(1.0-dx).*dy.*MTK(ym,xm);
            kt1(yn+1,xn)=kt1(yn+1,xn)+(1.0-dx).*dy.*MKT(MI(ym,xm));
            rhocp1(yn+1,xn)=rhocp1(yn+1,xn)+(1.0-dx).*dy.*MRHO(MI(ym,xm))*MCP(MI(ym,xm));
            wtnodes(yn+1,xn)=wtnodes(yn+1,xn)+(1.0-dx).*dy;
            rho1(yn,xn+1)=rho1(yn,xn+1)+dx.*(1.0-dy).*MRHO(MI(ym,xm));
            tk1(yn,xn+1)=tk1(yn,xn+1)+dx.*(1.0-dy).*MTK(ym,xm);
            kt1(yn,xn+1)=kt1(yn,xn+1)+dx.*(1.0-dy).*MKT(MI(ym,xm));
            rhocp1(yn,xn+1)=rhocp1(yn,xn+1)+dx.*(1.0-dy).*MRHO(MI(ym,xm))*MCP(MI(ym,xm));
            wtnodes(yn,xn+1)=wtnodes(yn,xn+1)+dx.*(1.0-dy);
            rho1(yn+1,xn+1)=rho1(yn+1,xn+1)+dx.*dy.*MRHO(MI(ym,xm));
            tk1(yn+1,xn+1)=tk1(yn+1,xn+1)+dx.*dy.*MTK(ym,xm);
            kt1(yn+1,xn+1)=kt1(yn+1,xn+1)+dx.*dy.*MKT(MI(ym,xm));
            rhocp1(yn+1,xn+1)=rhocp1(yn+1,xn+1)+dx.*dy.*MRHO(MI(ym,xm))*MCP(MI(ym,xm));
            wtnodes(yn+1,xn+1)=wtnodes(yn+1,xn+1)+dx.*dy;

            % Computing Marker Viscosity
            % 0 = Constant Viscosity
            METACUR=META(MI(ym,xm),1);
            % 1 = Isothermal power-law 2*EPSILONii=C1*SIGMAii^n
            % ETA=(2/C1)^(1/n)/2/EPSii^(1-1/n)
            if (MLAW(MI(ym,xm))==1)
                % ETA=(2/C1)^(1/n)/2/EPSii^(1-1/n)
                C1=META(MI(ym,xm),1);
                n=META(MI(ym,xm),2);
                meiicur=MEII(ym,xm);
                % Check second invariant of strain rate
                if (meiicur<1e-20) 
                    meiicur=1e-20;
                end
                METACUR=C1^(-1/n)*(2*meiicur)^(1/n-1);
            end
            % 2 = temperature dependent ETA=A*exp(Ea/RT0(1-(T-T0)/T0))
            if (MLAW(MI(ym,xm))==2)
                A=META(MI(ym,xm),1);
                Ea=META(MI(ym,xm),2);
                T0=META(MI(ym,xm),3);
                T=MTK(ym,xm);
                METACUR=A*exp(Ea/RGAS/T0*(1-(T-T0)/T0));
            end

            % Add shear viscosity etas() and rock type typ() to 4 surrounding nodes
            % only using markers located at <=0.5 gridstep distances from nodes
            if(dx<=0.5 && dy<=0.5)
                etas1(yn,xn)=etas1(yn,xn)+(1.0-dx).*(1.0-dy).*METACUR;
                typ1(yn,xn)=typ1(yn,xn)+(1.0-dx).*(1.0-dy).*MI(ym,xm);
                wtetas(yn,xn)=wtetas(yn,xn)+(1.0-dx).*(1.0-dy);
            end
            if(dx<=0.5 && dy>=0.5)
                etas1(yn+1,xn)=etas1(yn+1,xn)+(1.0-dx).*dy.*METACUR;
                typ1(yn+1,xn)=typ1(yn+1,xn)+(1.0-dx).*dy.*MI(ym,xm);
                wtetas(yn+1,xn)=wtetas(yn+1,xn)+(1.0-dx).*dy;
            end
            if(dx>=0.5 && dy<=0.5)
                etas1(yn,xn+1)=etas1(yn,xn+1)+dx.*(1.0-dy).*METACUR;
                typ1(yn,xn+1)=typ1(yn,xn+1)+dx.*(1.0-dy).*MI(ym,xm);
                wtetas(yn,xn+1)=wtetas(yn,xn+1)+dx.*(1.0-dy);
            end
            if(dx>=0.5 && dy>=0.5)
                etas1(yn+1,xn+1)=etas1(yn+1,xn+1)+dx.*dy.*METACUR;
                typ1(yn+1,xn+1)=typ1(yn+1,xn+1)+dx.*dy.*MI(ym,xm);
                wtetas(yn+1,xn+1)=wtetas(yn+1,xn+1)+dx.*dy;
            end

            % Add normal viscosity etan() to the center of current cell
            etan1(yn,xn)=etan1(yn,xn)+(1.0-abs(0.5-dx)).*(1.0-abs(0.5-dy)).*METACUR;
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
                tk1(i,j)=tk1(i,j)./wtnodes(i,j);
                kt1(i,j)=kt1(i,j)./wtnodes(i,j);
                rhocp1(i,j)=rhocp1(i,j)./wtnodes(i,j);
            else
                % If no new value is interpolated from markers old value is used
                rho1(i,j)=rho0(i,j);
                tk1(i,j)=tk0(i,j);
                kt1(i,j)=kt0(i,j);
                rhocp1(i,j)=rhocp0(i,j);
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
    
    % Applying thermal boundary conditions for interpolated temperature
    % Upper, Lower boundaries
    for j=2:1:xnum-1
        % Upper boundary
        tk1(1,j)=btopt(j,1)+btopt(j,2)*tk1(2,j);
        % Lower boundary
        tk1(ynum,j)=bbottomt(j,1)+bbottomt(j,2)*tk1(ynum-1,j);
    end
    % Left, Right boundaries: constant temperature
    for i=1:1:ynum
        % Left boundary
        tk1(i,1)=bleftt(i,1)+bleftt(i,2)*tk1(i,2);
        % Right boundary
        tk1(i,xnum)=brightt(i,1)+brightt(i,2)*tk1(i,xnum-1);
    end


    % Computing right part of mechanical equations
    % x-Stokes
    RX1=zeros(ynum+1,xnum);
    % y-Stokes
    RY1=zeros(ynum,xnum+1);
    % continuity
    RC1=zeros(ynum-1,xnum-1);
    % Grid points cycle
    for i=1:1:ynum;
        for j=1:1:xnum;
            % Right part of x-Stokes Equation
            if(j>1 && i>1 && j<xnum)
                RX1(i,j)=-gx*(rho1(i,j)+rho1(i-1,j))/2;
            end
            % Right part of y-Stokes Equation
            if(j>1 && i>1 && i<ynum)
                RY1(i,j)=-gy*(rho1(i,j)+rho1(i,j-1))/2;
            end
        end
    end


    % Solving of Stokes and Continuity equations on nodes
    % and computing residuals
    % by calling function Stokes_Continuity_solver_ghost()
    [vx1,resx1,vy1,resy1,pr1,resc1]=Stokes_Continuity_solver_Couette(prfirst,etas1,etan1,xnum,ynum,xstp,ystp,RX1,RY1,RC1,bleft,bright,btop,bbottom);

    % Computing EPS'xx=-EPS'yy, EPSxy=EPSyx deviatoric strain rate tensor components from vx, vy
    % Grid points cycle
    for i=1:1:ynum;
        for j=1:1:xnum;
            % EPS'xx=-EPS'yy=1/2(dvx/dx-dvy/dy)
            if(i<ynum & j<xnum)
                exx(i,j)=0.5*((vx1(i+1,j+1)-vx1(i+1,j))/xstp-(vy1(i+1,j+1)-vy1(i,j+1))/ystp);
            end
            % EPSxy=EPSyx=1/2(dvx/dy+dvy/dx)
            exy(i,j)=0.5*((vx1(i+1,j)-vx1(i,j))/ystp+(vy1(i,j+1)-vy1(i,j))/xstp);
        end
    end

    % Defining displacement timestep
    timestep=timemax
    % Check maximal velocity
    vxmax=max(abs(max(max(vx1))),abs(min(min(vx1))));
    vymax=max(abs(max(max(vy1))),abs(min(min(vy1))));
    % Check marker displacement step
    if (vxmax>0)
        if (timestep>markmax*xstp/vxmax);
            timestep=markmax*xstp/vxmax;
        end
    end
    if (vymax>0)
        if (timestep>markmax*ystp/vymax);
            timestep=markmax*ystp/vymax;
        end
    end
    timestep=timestep
    
    % Solving Temperature equation
    if (timestep>0 && tempmax>0)
        
        % Computing right part of temperature equation
        RT1=zeros(ynum,xnum);
        % Computing shear heating for Temperature nodes
        % Hs=4*etan*EPS'xx^2+4*etas*EPSxy^2
        % Grid points cycle
        for i=2:1:ynum-1;
            for j=2:1:xnum-1;
                % Adding 4*etas*EPSxy^2
                RT1(i,j)=RT1(i,j)+4*etas1(i,j)*exy(i,j)^2;
                % Computing and adding 4*etan*EPS'xx^2
                RT1(i,j)=RT1(i,j)+etan1(i-1,j-1)*exx(i-1,j-1)^2+etan1(i,j-1)*exx(i,j-1)^2+etan1(i-1,j)*exx(i-1,j)^2+etan1(i,j)*exx(i,j)^2;
            end
        end

        % Solving Temperature equation with displacement timestep
        [tk2,rest]=Temperature_solver(timestep,xnum,ynum,xstp,ystp,kt1,rhocp1,tk1,RT1,bleftt,brightt,btopt,bbottomt);
        % Computing temperature changes
        dtk1=tk2-tk1;
        % Checking temperature changes
        dtkmax=max(max(abs(dtk1)))
        % Repeating temperature solution if temperature changes are too big
        if(dtkmax>tempmax)
            % Computing reduced timestep
            timestep=timestep*tempmax/dtkmax;
            % Solving Temperature equation with reduced timestep
            [tk2,rest]=Temperature_solver(timestep,xnum,ynum,xstp,ystp,kt1,rhocp1,tk1,RT1,bleftt,brightt,btopt,bbottomt);
            % Computing temperature changes
            dtk1=tk2-tk1;
        end

        % Computing subgrid diffusion for markers
        if (dsubgridt>0)
            % Clear subgrid temperature changes for nodes
            dtkn=zeros(ynum,xnum);
            % Clear wights for basic nodes
            wtnodes=zeros(ynum,xnum);
            % Marker cycle
            for xm = 1:1:mxnum
                for ym = 1:1:mynum  
                    %  yn    T(yn,xn)--------------------T(yn,xn+1)
                    %           ?           ^                  ?
                    %           ?           ?                  ?
                    %           ?          dy                  ?
                    %           ?           ?                  ?
                    %           ?           v                  ?
                    %           ?<----dx--->o Mrho(ym,xm)       ?
                    %           ?                              ?
                    %           ?                              ?
                    %  yn+1  T(yn+1,xn)-------------------V(yn+1,xn+1)
                    %
                    %
                    % Interpolating temperature changes from basic nodes
                    % !!! SUBTRACT 0.5 since int16(0.5)=1
                    xn=double(int16(MX(ym,xm)./xstp-0.5))+1;
                    yn=double(int16(MY(ym,xm)./ystp-0.5))+1;
                    % Check horizontal index for upper left T-node 
                    % It must be between 1 and xnum-1 (see picture for grid)
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
                    % Define and check normalized distances from marker to the upper left T-node;
                    dx=MX(ym,xm)./xstp-xn+1;
                    dy=MY(ym,xm)./ystp-yn+1;
                    % Interpolate old nodal temperature for the marker
                    tkm=0;
                    tkm=tkm+(1.0-dx).*(1.0-dy).*tk1(yn,xn);
                    tkm=tkm+(1.0-dx).*dy.*tk1(yn+1,xn);
                    tkm=tkm+dx.*(1.0-dy).*tk1(yn,xn+1);
                    tkm=tkm+dx.*dy.*tk1(yn+1,xn+1);
                    % Calculate Nodal-Marker subgrid temperature difference
                    dtkm=tkm-MTK(ym,xm);
                    % Compute nodal k and RHO*Cp for the marker 
                    % k
                    ktm=0;
                    ktm=ktm+(1.0-dx).*(1.0-dy).*kt1(yn,xn);
                    ktm=ktm+(1.0-dx).*dy.*kt1(yn+1,xn);
                    ktm=ktm+dx.*(1.0-dy).*kt1(yn,xn+1);
                    ktm=ktm+dx.*dy.*kt1(yn+1,xn+1);
                    % RHO*Cp
                    rhocpm=0;
                    rhocpm=rhocpm+(1.0-dx).*(1.0-dy).*rhocp1(yn,xn);
                    rhocpm=rhocpm+(1.0-dx).*dy.*rhocp1(yn+1,xn);
                    rhocpm=rhocpm+dx.*(1.0-dy).*rhocp1(yn,xn+1);
                    rhocpm=rhocpm+dx.*dy.*rhocp1(yn+1,xn+1);
                    
                    % Compute local thermal diffusion timescale for the marker
                    tdm=rhocpm/ktm/(2/xstp^2+2/ystp^2);
                    
                    % Computing subgrid diffusion
                    sdif=-dsubgridt*timestep/tdm;
                    if(sdif<-30) 
                        sdif=-30;
                    end
                    dtkm=dtkm*(1-exp(sdif));    
                    
                    % Correcting old temperature for the marker
                    MTK(ym,xm)=MTK(ym,xm)+dtkm;
                    
                    % Interpolating subgrid temperature changes to 4
                    dtkn(yn,xn)=dtkn(yn,xn)+(1.0-dx).*(1.0-dy).*dtkm;
                    wtnodes(yn,xn)=wtnodes(yn,xn)+(1.0-dx).*(1.0-dy);
                    dtkn(yn+1,xn)=dtkn(yn+1,xn)+(1.0-dx).*dy.*dtkm;
                    wtnodes(yn+1,xn)=wtnodes(yn+1,xn)+(1.0-dx).*dy;
                    dtkn(yn,xn+1)=dtkn(yn,xn+1)+dx.*(1.0-dy).*dtkm;
                    wtnodes(yn,xn+1)=wtnodes(yn,xn+1)+dx.*(1.0-dy);
                    dtkn(yn+1,xn+1)=dtkn(yn+1,xn+1)+dx.*dy.*dtkm;
                    wtnodes(yn+1,xn+1)=wtnodes(yn+1,xn+1)+dx.*dy;
                    
                end
            end

            % Computing subgrid diffusion for nodes
            for i=1:1:ynum;
                for j=1:1:xnum;
                    % Density
                    if (wtnodes(i,j)~=0)
                        % Compute new value interpolated from markers
                        dtkn(i,j)=dtkn(i,j)./wtnodes(i,j);
                    end
                end
            end
            
            % Subtracting subgrid diffusion part from nodal temperature changes
            dtk1=dtk1-dtkn;
            
        end
        
        % Updating temperature for markers
        for xm = 1:1:mxnum
            for ym = 1:1:mynum  

                %  yn    T(yn,xn)--------------------T(yn,xn+1)
                %           ?           ^                  ?
                %           ?           ?                  ?
                %           ?          dy                  ?
                %           ?           ?                  ?
                %           ?           v                  ?
                %           ?<----dx--->o Mrho(ym,xm)       ?
                %           ?                              ?
                %           ?                              ?
                %  yn+1  T(yn+1,xn)-------------------V(yn+1,xn+1)
                %
                %
                % Interpolating temperature changes from basic nodes
                % !!! SUBTRACT 0.5 since int16(0.5)=1
                xn=double(int16(MX(ym,xm)./xstp-0.5))+1;
                yn=double(int16(MY(ym,xm)./ystp-0.5))+1;
                % Check horizontal index for upper left T-node 
                % It must be between 1 and xnum-1 (see picture for grid)
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
                % Define and check normalized distances from marker to the upper left T-node;
                dx=MX(ym,xm)./xstp-xn+1;
                dy=MY(ym,xm)./ystp-yn+1;
                % Calculate Marker temperature change from four surrounding nodes
                dtkm=0;
                dtkm=dtkm+(1.0-dx).*(1.0-dy).*dtk1(yn,xn);
                dtkm=dtkm+(1.0-dx).*dy.*dtk1(yn+1,xn);
                dtkm=dtkm+dx.*(1.0-dy).*dtk1(yn,xn+1);
                dtkm=dtkm+dx.*dy.*dtk1(yn+1,xn+1);
                %
                %Computing new temperature for the marker
                MTK(ym,xm)=MTK(ym,xm)+dtkm;
            end
        end

    end
    
    

    % Computing strain rate invariant for markers
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


            % Calculating Second invariant of strain rate for marker
            %
            % Interpolating squares of EPS'xx=-EPS'yy from cell centers
            % EPS'xx-nodes are displaced rightward and downward for 1/2 of gridsteps
            % !!! SUBTRACT 0.5 since int16(0.5)=1
            xn=double(int16((MX(ym,xm)-xstp/2.0)./xstp-0.5))+1;
            yn=double(int16((MY(ym,xm)-ystp/2.0)./ystp-0.5))+1;
            % Check horizontal index for upper left VY-node 
            % It must be between 1 and xnum-2 (see picture for staggered grid)
            if (xn<1)
                xn=1;
            end
            if (xn>xnum-2)
                xn=xnum-2;
            end
            if (yn<1)
                yn=1;
            end
            if (yn>ynum-2)
                yn=ynum-2;
            end
            % Define and check normalized distances from marker to the upper left EPSVX-node;
            dx=(MX(ym,xm)-xstp/2.0)./xstp-xn+1;
            dy=(MY(ym,xm)-ystp/2.0)./ystp-yn+1;
            % Calculate Marker EPS'xx^2 from four surrounding nodes
            exxm=0;
            exxm=exxm+(1.0-dx).*(1.0-dy).*exx(yn,xn)*exx(yn,xn);
            exxm=exxm+(1.0-dx).*dy.*exx(yn+1,xn)*exx(yn+1,xn);
            exxm=exxm+dx.*(1.0-dy).*exx(yn,xn+1)*exx(yn,xn+1);
            exxm=exxm+dx.*dy.*exx(yn+1,xn+1)*exx(yn+1,xn+1);
            %
            % Interpolating squares of EPSxy=EPSyx from basic nodes
            % !!! SUBTRACT 0.5 since int16(0.5)=1
            xn=double(int16(MX(ym,xm)./xstp-0.5))+1;
            yn=double(int16(MY(ym,xm)./ystp-0.5))+1;
            % Check horizontal index for upper left VY-node 
            % It must be between 1 and xnum-2 (see picture for staggered grid)
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
            % Define and check normalized distances from marker to the upper left VX-node;
            dx=MX(ym,xm)./xstp-xn+1;
            dy=MY(ym,xm)./ystp-yn+1;
            % Calculate Marker EPSxy^2 from four surrounding nodes
            exym=0;
            exym=exym+(1.0-dx).*(1.0-dy).*exy(yn,xn)*exy(yn,xn);
            exym=exym+(1.0-dx).*dy.*exy(yn+1,xn)*exy(yn+1,xn);
            exym=exym+dx.*(1.0-dy).*exy(yn,xn+1)*exy(yn,xn+1);
            exym=exym+dx.*dy.*exy(yn+1,xn+1)*exy(yn+1,xn+1);
            %
            %Computing second strain rate invariant for the marker
            MEII(ym,xm)=(exxm+exym)^0.5;
       end
    end


    % Moving Markers by velocity field
    if(markmove==1)
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

                % Define indexes for upper left node in the VX-cell where the marker is
                % VX-cells are displaced upward for 1/2 of vertical gridstep
                % !!! SUBTRACT 0.5 since int16(0.5)=1
                xn=double(int16(MX(ym,xm)./xstp-0.5))+1;
                yn=double(int16((MY(ym,xm)+ystp/2.0)./ystp-0.5))+1;
                % Check vertical index for upper left VX-node 
                % It must be between 1 and ynum-2 (see picture for staggered grid)
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
                % It must be between 1 and xnum-2 (see picture for staggered grid)
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

                % Recycling markers exiting the model
                % Top->Bottom
                if (MY(ym,xm)<0)
                    MY(ym,xm)=ysize+MY(ym,xm);
                    MTK(ym,xm)=MTK(ym,xm)-ttop+tbottom;
                end
                 % Bottom->Top
                if (MY(ym,xm)>ysize)
                    MY(ym,xm)=MY(ym,xm)-ysize;
                    MTK(ym,xm)=MTK(ym,xm)+ttop-tbottom;
                end
           end
        end
    end

    % Defining scale for Stokes residuals from y-Stokes equation
    g=10;
    % dSIGMAij/dj-dP/di=-RHO*gi=0  => Stokes scale=abs(RHO*gi)
    stokesscale= MRHO(1)*g;
    % Defining scale for Continuity residuals from y-Stokes equation
    % dvx/dx+dvy/dy=0 can be transformed to 2ETA(dvx/dx+dvy/dy)/dx=0 
    % which is similar to dSIGMAij/dj and has scale given above
    % therefore continuity scale = scale=abs(RHO*gi/ETA*dx)
    continscale= MRHO(1)*g/etas1(1,1)*xstp;


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
    surf(tk1);
    shading interp;
    light;
    lighting phong;
    axis tight;
    zlabel(['T K, Time = ',num2str(timesum*1e-6/(365.25*24*3600)),' Myr']);
    
    vy0=vy1*1E+10;
    
    timesum=timesum
    timestep=timestep
    % Advance in time
    timesum=timesum+timestep
    
    % Computing FI1
    T1=tk2(1,xnum)
    TET1=META(1,2)*(T1-META(1,3))/RGAS/META(1,3)^2
    FI1=exp(TET1)
    % Computing B parameter for given FI1 by bisection
    % Lower limit
    Bmin=1e-6;
    BRmin=Bmin^2/2*(1-((exp(Bmin)-1)/(exp(Bmin)+1))^2);
    FI1min=Bmin^2/2/BRmin;
    dFI1min=FI1min-FI1;
    % Upper limit
    Bmax=10;
    BRmax=Bmax^2/2*(1-((exp(Bmax)-1)/(exp(Bmax)+1))^2);
    FI1max=Bmax^2/2/BRmax;
    dFI1max=FI1max-FI1;
    % Bisection
    if (dFI1min<=0 && dFI1max>=0)
        dFI1cur=1e+6;
        while(abs(dFI1cur)>1e-6)
            Bcur=(Bmin+Bmax)/2;
            BRcur=Bcur^2/2*(1-((exp(Bcur)-1)/(exp(Bcur)+1))^2);
            FI1cur=Bcur^2/2/BRcur;
            dFI1cur=FI1cur-FI1;
            if(dFI1cur<=0)
                Bmin=Bcur;
            else
                Bmax=Bcur;
            end
        end
    end
    
    % Brinkman number based on temperature
    BRT(ntimestep)=BRcur;
    BRT1=BRT(ntimestep)

    % Brinkman number based on stress
    SIGyx=2*exy(2,2)*etas1(2,2)
    BRS(ntimestep)=(xsize*SIGyx)^2*META(1,2)/(MKT(1)*META(1,1)*RGAS*META(1,3)^2)/exp(META(1,2)/RGAS/META(1,3));
    BRS1=BRS(ntimestep)
    
   
    % Br difference
    BRD(ntimestep)=log10(abs(BRS(ntimestep)-BRT(ntimestep)));
    % Plotting Brinkman number
    subplot(2,3,4)
    plot(log10(BRS),'k');
    hold on;
    plot(log10(BRT),'r');
    hold on;
    plot(BRD,'b');
    ylabel(['Br convergence, Step = ',num2str(ntimestep)])
%     colorbar;
    
    

    pause(0.01)

end