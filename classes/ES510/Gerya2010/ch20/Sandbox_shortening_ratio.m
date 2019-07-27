% Solution of Stokes and continuity equations for viscoelastoplastic medium 
% with variable viscosity and variable shear modulus in 2D with direct solver
% by using external function Stokes_Continuity_solver_sandbox()
% Setup corresponds to viscoelastoplastic block deformation
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

% Temperature at the top, and bottom of the model, K
ttop=273;
tbottom=273;


% Acceleration of Gravity, m/s^2
gx=0;
gy=9.81;

% Horizontal extension strain rate, dvx/dx, 1/s
epsext=0e-14;

% Gas constant J/mol/K
RGAS=8.314;

% Model size, m
xsize=0.19;
ysize=0.06;

% Defining resolution
xnum=191;
ynum=61;

% Viscoelastic medium
% Rock density (kg/m3): RHO*exp(-ALP*(T-T0))
MRHO(1,1)=1;               % standard density, kg/m^3
MRHO(1,2)=0;             % thermal expansion, 1/K
MRHO(1,3)=ttop;   % T0
MRHO(2,1)=1560;               % standard density, kg/m^3
MRHO(2,2)=0;             % thermal expansion, 1/K
MRHO(2,3)=ttop;   % T0
MRHO(3,1)=1560;               % standard density, kg/m^3
MRHO(3,2)=0;             % thermal expansion, 1/K
MRHO(3,3)=ttop;   % T0% Rock flow law: 
MRHO(4,1)=1480;               % standard density, kg/m^3
MRHO(4,2)=0;             % thermal expansion, 1/K
MRHO(4,3)=ttop;   % T0
MRHO(5,1)=1520;               % standard density, kg/m^3
MRHO(5,2)=0;             % thermal expansion, 1/K
MRHO(5,3)=ttop;   % T0
% 0 = constant viscosity;
% 1 = isothermal power-law 2*EPSILONii=C1*SIGMAii^n
% 2 = temperature dependent ETA=A*exp(Ea/RT0(1-(T-T0)/T0))
% 3 = temperature & depth dependent ETA=rho*A*exp(-b*(T-ttop)/(tbottom-ttop) + c*Y/yzize)
MFLOW(1,1)=0;
MFLOW(1,2)=1e+2;     % ETA, Pa s
MFLOW(2,1)=0;
MFLOW(2,2)=1e+09;     % ETA, Pa s
MFLOW(3,1)=0;
MFLOW(3,2)=1e+09;     % ETA, Pa s
MFLOW(4,1)=0;
MFLOW(4,2)=1e+09;     % ETA, Pa s
MFLOW(5,1)=0;
MFLOW(5,2)=1e+12;     % ETA, Pa s
% Viscosity limits for the sand
etamin=1e+4;
etamax=1e+9;
% Rock shear modulus, Pa
MMU(1)=3e+6;
MMU(2)=3e+6;
MMU(3)=3e+6;
MMU(4)=3e+6;
MMU(5)=3e+16;
% Mohr-Coulomb plasticity: SIGMAyeild=C+sin(FI)*P
% C=C0, FI=FI0 for strain<=GAM0
% C=C0+(C1-C0)/(GAM1-GAM0)*(strain-GAM0), FI=FI0+(FI1-FI0)/(GAM1-GAM0)*(strain-GAM0) for GAM0<strain<GAM1
% C=C1, FI=FI1 for strain>=GAM0
% AIR
MPL(1,1)=1e+38; % C0, Pa
MPL(1,2)=1e+38; % C1, Pa
MPL(1,3)=0;     % sin(FI0)
MPL(1,4)=0;     % sin(FI1)
MPL(1,5)=0;     % GAM0
MPL(1,6)=1;     % GAM1
% Sand1
MPL(2,1)=10; % C0, Pa
MPL(2,2)=10; % C1, Pa
MPL(2,3)=sin(36/180*pi);     % sin(FI0)
MPL(2,4)=sin(31/180*pi);     % sin(FI1)
MPL(2,5)=0;     % GAM0
MPL(2,6)=1;     % GAM1
% Sand2
MPL(3,1)=10; % C0, Pa
MPL(3,2)=10; % C1, Pa
MPL(3,3)=sin(36/180*pi);     % sin(FI0)
MPL(3,4)=sin(31/180*pi);     % sin(FI1)
MPL(3,5)=0;     % GAM0
MPL(3,6)=1;     % GAM1
% Microbeads
MPL(4,1)=10; % C0, Pa
MPL(4,2)=10; % C1, Pa
MPL(4,3)=sin(22/180*pi);     % sin(FI0)
MPL(4,4)=sin(20/180*pi);     % sin(FI1)
MPL(4,5)=0;     % GAM0
MPL(4,6)=1;     % GAM1
% Moving wall
MPL(5,1)=1e+38; % C0, Pa
MPL(5,2)=1e+38; % C1, Pa
MPL(5,3)=0;     % sin(FI0)
MPL(5,4)=0;     % sin(FI1)
MPL(5,5)=0;     % GAM0
MPL(5,6)=1;     % GAM1
% Rock heat capacity, J/K/kg
MCP(1)=1000;
MCP(2)=1000;
MCP(3)=1000;
MCP(4)=1000;
MCP(5)=1000;
% Rock thermal conductivity, W/m/K
MKT(1,1)=3; 
MKT(2,1)=3; 
MKT(3,1)=3; 
MKT(4,1)=3; 
MKT(5,1)=3; 
% Viscoelastic timestep, s
timemax=6e-2;
% Maximal marker displacement step, number of gridsteps
markmax=0.1;
% Moving Markers: 
% 0 = not moving at all
% 1 = simple advection
% 4 = 4-th order in space  Runge-Kutta
markmove=4;
% Velocity calculation
% 0 = by Solving momentum and continuity equations
% 1 = solid body rotation
movemod=0;
% Maximal temperature change, allowed for one timestep, K
tempmax=0;
% Amount of timesteps
stepmax=10000;

% Numerical Subgrid stress diffusion coefficient
dsubgrids=1;
% Numerical Subgrid temperature diffusion coefficient
dsubgridt=1;
% Shear heating on(1)/off(0)
frictyn=0;
% Weight for old visco-plastic viscosity
etawt=0.0; 

% Pressure boundary conditions
% prfirst(1) = boundary condition mode:
% 0 - pressure in one cell definition
% 1 - pressure at the top and in the bottom of the channel
prfirst(1)=0;
% prfirst(2) = boundary condition value
prfirst(2)=20;


% Velocity Boundary condition specified by bleft,bright,btop,bbot 
% are implemented from ghost nodes 
% directly into Stokes and continuity equations

% Upper, Lower boundaries: Free slip + Prescribed inward velocity (vertical shortening)
for j=1:1:xnum+1
    % Upper boundary: Free slip
    % vx(1,j)=btop(j,1)+vx(2,j)*btop(j,2)
    btop(j,1)=0;
    btop(j,2)=1;
    % vy(1,j)=btop(j,3)+vy(2,j)*btop(j,4)
    btop(j,3)=epsext*ysize/2;
    btop(j,4)=0;
    % Lower boundary: No Slip  
    % vx(ynum+1,j)=bbottom(j,1)+vx(ynum,j)*bbottom(j,2)
    bbottom(j,1)=0;
    bbottom(j,2)=-1;
    % vy(ynum,j)=bbottom(j,3)+vy(ynum-1,j)*bbottom(j,4)
    bbottom(j,3)=-epsext*ysize/2;
    bbottom(j,4)=0;
end

% Left, Right boundaries: + Prescribed outward velocity (horizontal extension)
for i=1:1:ynum+1
    % Left boundary: No slip   
    % vx(i,1)=bleft(i,1)+vx(i,2)*bleft(i,2)
    bleft(i,1)=-epsext*xsize/2;
    bleft(i,2)=0;
    % vy(i,1)=bleft(i,3)+vy(i,2)*bleft(i,42)
    bleft(i,3)=0;
    bleft(i,4)=-1;
    % Right boundary: Free slip 
    % vx(i,xnum)=bright(i,1)+vx(i,xnum-1)*bbright(i,2)
    bright(i,1)=epsext*xsize/2;
    bright(i,2)=0;
    % vy(i,xnum+1)=bright(i,3)+vx(i,xnum)*bbright(i,4)
    bright(i,3)=0;
    bright(i,4)=1;
end

% Internal boundary condition: prescribed velocity of mobile wall
bintern(1)=xnum-8;  % Horizontal position of vx nodes with prescrbed velocity 
bintern(2)=6;       % Min vertical position
bintern(3)=ynum-3;  % Max vertical position
bintern(4)=-2.5/100/3600; % Prescribed shortening velocity, m/s (2.5 cm/hour)
bintern(5)=xnum-8;  % Horizontal position of vy nodes with prescrbed velocity 
bintern(6)=ynum-3;  % Min vertical position
bintern(7)=ynum-3;  % Max vertical position
bintern(8)=0; % No vertical velocity in the wall bottom

% Defining average gridsteps
xstp=xsize./(xnum-1);
ystp=ysize./(ynum-1);

% Defining gridline positions for regular basic grid
% In horizontal direction
gridx=zeros(xnum,1);
for i=2:1:xnum
    gridx(i)=gridx(i-1)+xstp;
end
% In vertical direction
gridy=zeros(ynum,1);
for i=2:1:ynum
    gridy(i)=gridy(i-1)+ystp;
end


% Thermal boundary conditions
% Upper, Lower boundaries: constant temperature
for j=1:1:xnum
    % Upper boundary
    % tk(1,j)=btopt(j,1)+tk(2,j)*btop(j,2)
    btopt(j,1)=ttop;
    btopt(j,2)=0;
    % Lower boundary
    % tk(ynum,j)=bbottomt(j,1)+tk(ynum-1,j)*bbottomt(j,2)
    bbottomt(j,1)=tbottom;
    bbottomt(j,2)=0;
end
% Left, Right boundaries: symmetry
for i=1:1:ynum
    % Left boundary
    % tk(i,1)=bleftt(i,1)+bleftt(i,2)*tk(i,2);
    bleftt(i,1)=0;
    bleftt(i,2)=1;
    % Right boundary
    % tk(i,xnum)=brightt(i,1)+brightt(i,2)*tk(i,xnum-1);
    brightt(i,1)=0;
    brightt(i,2)=1;
end


% Defining number of markers and steps between them in the horizontal and vertical direction
mxnum=760; %total number of markers in horizontal direction
mynum=240; %total number of markers in vertical direction
mxstep=xsize/mxnum; %step between markers in horizontal direction   
mystep=ysize/mynum; %step between markers in vertical direction

% Creating markers arrays
MX=zeros(mynum,mxnum);   % X coordinate, m
MY=zeros(mynum,mxnum);   % Y coordinate, m
MTK=zeros(mynum,mxnum);  % Temperature, K
MI=zeros(mynum,mxnum);   % Type
MXN=zeros(mynum,mxnum);  % Horizontal index
MYN=zeros(mynum,mxnum);  % Vertical index
MSXX=zeros(mynum,mxnum);  % SIGMAxx - deviatoric normal stress, Pa
MSXY=zeros(mynum,mxnum);  % SIGMAyy - shear stress, Pa
META=zeros(mynum,mxnum);  % viscosity, Pa s
MEXX=zeros(mynum,mxnum);  % EPSILONxx - normal strain rate, 1/s
MEXY=zeros(mynum,mxnum);  % EPSILONyy - shear strain rate, 1/s
MPR=zeros(mynum,mxnum);   % Pressure, Pa
MGII=zeros(mynum,mxnum);  % Accumulated strain
MRAT=ones(mynum,mxnum);  % EiiMarker/EiiGrid Ratio

% Defining intial position of markers
% Defining lithological structure of the model
for xm = 1:1:mxnum
    for ym = 1:1:mynum
        % Coordinates with small random displacement
        MX(ym,xm)=xm*mxstep-mxstep/2+(rand-0.5)*mxstep;
        MY(ym,xm)=ym*mystep-mystep/2+(rand-0.5)*mystep;
        % Material Type
        MI(ym,xm)=2;
        % Layered structure in the sand
        m2=double(int16(MY(ym,xm)/(ysize/24)-0.5))+1;
        m2=m2-double(int16(m2/2-0.5))*2;
        if(m2==0)
            MI(ym,xm)=3;
        end
        % Microbeads
        if (MY(ym,xm)<ysize-0.005 && MY(ym,xm)>ysize-0.010)
            MI(ym,xm)=4;
        end
        % Mobile Wall
        if (MX(ym,xm)>xsize-0.01 && MX(ym,xm)<xsize-0.005 && MY(ym,xm)>0.005 && MY(ym,xm)<ysize-0.002)
            MI(ym,xm)=5;
        end
        % Air
        if (MX(ym,xm)<=xsize-0.11 && MY(ym,xm)<=ysize-0.035)
            MI(ym,xm)=1;
        end
        if (MX(ym,xm)>=xsize-0.005 && MY(ym,xm)<=ysize-0.002)
            MI(ym,xm)=1;
        end
        if (MY(ym,xm)<=0.005)
            MI(ym,xm)=1;
        end
        if (MX(ym,xm)>xsize-0.11 && MX(ym,xm)<=xsize-0.01 && MY(ym,xm)<=ysize-0.035-(MX(ym,xm)-(xsize-0.11))*tan(10/180*pi))
            MI(ym,xm)=1;
        end
        % Initial temperature profile
        ywt=MY(ym,xm)/ysize;
        MTK(ym,xm)=ttop+(tbottom-ttop)*ywt;
    end
end


% Density, viscosity, shear modulus, temperature, thermal conductivity, RHO*Cp arrays
etas1 = zeros(ynum,xnum);       % Viscosity for shear stress
etan1 = zeros(ynum-1,xnum-1);   % Viscosity for normal stress
mus1 = zeros(ynum,xnum);        % Shear modulus for shear stress
mun1 = zeros(ynum-1,xnum-1);    % Shear modulus for normal stress
sxy1 = zeros(ynum,xnum);        % Shear stress
sxx1 = zeros(ynum-1,xnum-1);    % Normal stress
rho1 = zeros(ynum,xnum);        % Density
tk1 = zeros(ynum,xnum);         % Old temperature
tk2=tk1;                        % New temperature
rhocp1 = zeros(ynum,xnum);      % RHO*Cp (for temperature equation)
kt1 = zeros(ynum,xnum);         % Thermal conductivity

% Creating matrix for Stokes and continuity equations
L=sparse((xnum-1)*(ynum-1)*3,(xnum-1)*(ynum-1)*3);
R=zeros((xnum-1)*(ynum-1)*3,1);
% Creating matrix for temperature equation
LT=sparse(xnum*ynum,xnum*ynum);
RT=zeros(xnum*ynum,1);

% Initial time, s
timesum=0;
timestepd=timemax
ntimestep=1;
% Main Time cycle
for ntimestep=ntimestep:1:stepmax
    
    % Defining viscoelastic timestep
    timestep=timemax % viscoelastic timestep
    % Plastic yeilding mark
    plastyn=0;

    % Backup transport properties arrays
    etas0 = etas1;
    etan0 = etan1;
    mus0 = mus1;
    mun0 = mun1;
    sxy0 = sxy1;
    sxx0 = sxx1;
    rho0 = rho1;
    tk0=tk2;
    rhocp0=rhocp1;
    kt0=kt1;
    % Clear transport properties arrays
    etas1 = zeros(ynum,xnum);
    etan1 = zeros(ynum-1,xnum-1);
    mus1 = zeros(ynum,xnum);
    mun1 = zeros(ynum-1,xnum-1);
    sxy1 = zeros(ynum,xnum);
    sxx1 = zeros(ynum-1,xnum-1);
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

    % Computing grid steps for basic nodes
    xstp1=zeros(xnum-1,1);
    ystp1=zeros(ynum-1,1);
    for i=1:1:xnum-1
        xstp1(i)=gridx(i+1)-gridx(i);
    end
    for i=1:1:ynum-1
        ystp1(i)=gridy(i+1)-gridy(i);
    end  

    % Computing grids and grid steps for Vx, Vy nodes
    % Horizontal (for Vy)
    gridcx=zeros(xnum+1,1);
    xstpc1=zeros(xnum,1);
    % Vertical (for Vx)
    gridcy=zeros(ynum+1,1);
    ystpc1=zeros(ynum,1);
    % First and last nodes and steps (for external nodes)
    % Horizontal (for Vy)
    gridcx(1)=-xstp1(1)/2;
    xstpc1(1)=xstp1(1);
    gridcx(xnum+1)=xsize+xstp1(xnum-1)/2;
    xstpc1(xnum)=xstp1(xnum-1);
    % Vertical (for Vx)
    gridcy(1)=-ystp1(1)/2;
    ystpc1(1)=ystp1(1);
    gridcy(ynum+1)=ysize+ystp1(ynum-1)/2;
    ystpc1(ynum)=ystp1(ynum-1);
    % Internal nodes
    for i=2:1:xnum
        gridcx(i)=(gridx(i)+gridx(i-1))/2;
    end
    for i=2:1:ynum
        gridcy(i)=(gridy(i)+gridy(i-1))/2;
    end   
    % Internal grid steps
    for i=2:1:xnum-1
        xstpc1(i)=(gridx(i+1)-gridx(i-1))/2;
    end
    for i=2:1:ynum-1
        ystpc1(i)=(gridy(i+1)-gridy(i-1))/2;
    end   

    % Interpolating parameters from markers to nodes
    for xm = 1:1:mxnum
        for ym = 1:1:mynum 
            if(MI(ym,xm)==5)
                    MGII(ym,xm)=0;
            end
            % Check markers inside the grid
            if (MX(ym,xm)>=0 && MX(ym,xm)<=xsize && MY(ym,xm)>=0 && MY(ym,xm)<=ysize) 
                
                % Changing marker type to friction layer
                if((MI(ym,xm)==1 || MI(ym,xm)==5) && MY(ym,xm)>=ysize-0.002)
                    MI(ym,xm)=3;
                end
                % Changing marker type to growing mobile wall
                if(MI(ym,xm)~=5 && MX(ym,xm)>xsize-0.01 && MX(ym,xm)<xsize-0.005 && MY(ym,xm)<ysize-0.002 && MY(ym,xm)>0.005)
                    MI(ym,xm)=5;
                    MGII(ym,xm)=0;
                end

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
                % using bisection
                % Find horizontal index
                xnmin=1;
                xnmax=xnum;
                while ((xnmax-xnmin)>1)
                    % !!! SUBTRACT 0.5 since int16(0.5)=1
                    xn=double(int16((xnmax+xnmin)./2-0.5));
                    if(gridx(xn)>MX(ym,xm))
                        xnmax=xn;
                    else
                        xnmin=xn;
                    end
                end
                xn=xnmin;
                % Check horizontal index
                if (xn<1)
                    xn=1;
                end
                if (xn>xnum-1)
                    xn=xnum-1;
                end
                % Save horizontal index
                MXN(ym,xm)=xn;

                % Find vertical index
                ynmin=1;
                ynmax=ynum;
                while ((ynmax-ynmin)>1)
                    % !!! SUBTRACT 0.5 since int16(0.5)=1
                    yn=double(int16((ynmax+ynmin)./2-0.5));
                    if(gridy(yn)>MY(ym,xm))
                        ynmax=yn;
                    else
                        ynmin=yn;
                    end
                end
                yn=ynmin;
                % Check vertical index
                if (yn<1)
                    yn=1;
                end
                if (yn>ynum-1)
                    yn=ynum-1;
                end
                % Save Vertical index
                MYN(ym,xm)=yn;

                % Define normalized distances from marker to the upper left node;
                dx=(MX(ym,xm)-gridx(xn))./xstp1(xn);
                dy=(MY(ym,xm)-gridy(yn))./ystp1(yn);

                % Compute marker weight koefficient from cell dimensions
                % Number of markers in a cell is in invert proportion to the cell volume
                mwt=1;%/xstp1(xn)/ystp1(yn);

                % Compute density from marker temperature
                MRHOCUR=MRHO(MI(ym,xm),1)*(1-MRHO(MI(ym,xm),2)*(MTK(ym,xm)-ttop));

                % Compute rho*Cp for marker 
                %MRHOCPCUR=MRHOCUR*MCP(MI(ym,xm));
                MRHOCPCUR=MRHO(MI(ym,xm),1)*MCP(MI(ym,xm));

                % Compute thermal conductivity from marker temperature
                % Rock thermal conductivity, k=k0*rho*cp, W/m/K
                MKTCUR=MKT(MI(ym,xm),1)*MRHOCPCUR;

                % Add properties to 4 surrounding nodes
                rho1(yn,xn)=rho1(yn,xn)+(1.0-dx).*(1.0-dy).*MRHOCUR*mwt;
                tk1(yn,xn)=tk1(yn,xn)+(1.0-dx).*(1.0-dy).*MTK(ym,xm)*mwt;
                kt1(yn,xn)=kt1(yn,xn)+(1.0-dx).*(1.0-dy).*MKTCUR*mwt;
                rhocp1(yn,xn)=rhocp1(yn,xn)+(1.0-dx).*(1.0-dy).*MRHOCPCUR*mwt;
                wtnodes(yn,xn)=wtnodes(yn,xn)+(1.0-dx).*(1.0-dy)*mwt;

                rho1(yn+1,xn)=rho1(yn+1,xn)+(1.0-dx).*dy.*MRHOCUR*mwt;
                tk1(yn+1,xn)=tk1(yn+1,xn)+(1.0-dx).*dy.*MTK(ym,xm)*mwt;
                kt1(yn+1,xn)=kt1(yn+1,xn)+(1.0-dx).*dy.*MKTCUR*mwt;
                rhocp1(yn+1,xn)=rhocp1(yn+1,xn)+(1.0-dx).*dy.*MRHOCPCUR*mwt;
                wtnodes(yn+1,xn)=wtnodes(yn+1,xn)+(1.0-dx).*dy*mwt;

                rho1(yn,xn+1)=rho1(yn,xn+1)+dx.*(1.0-dy).*MRHOCUR*mwt;
                tk1(yn,xn+1)=tk1(yn,xn+1)+dx.*(1.0-dy).*MTK(ym,xm)*mwt;
                kt1(yn,xn+1)=kt1(yn,xn+1)+dx.*(1.0-dy).*MKTCUR*mwt;
                rhocp1(yn,xn+1)=rhocp1(yn,xn+1)+dx.*(1.0-dy).*MRHOCPCUR*mwt;
                wtnodes(yn,xn+1)=wtnodes(yn,xn+1)+dx.*(1.0-dy)*mwt;

                rho1(yn+1,xn+1)=rho1(yn+1,xn+1)+dx.*dy.*MRHOCUR*mwt;
                tk1(yn+1,xn+1)=tk1(yn+1,xn+1)+dx.*dy.*MTK(ym,xm)*mwt;
                kt1(yn+1,xn+1)=kt1(yn+1,xn+1)+dx.*dy.*MKTCUR*mwt;
                rhocp1(yn+1,xn+1)=rhocp1(yn+1,xn+1)+dx.*dy.*MRHOCPCUR*mwt;
                wtnodes(yn+1,xn+1)=wtnodes(yn+1,xn+1)+dx.*dy*mwt;

                % Computing Marker Viscosity
                METACUR=MFLOW(MI(ym,xm),2);

                % Check if any plastic yeiding condition is present 
                % Checking for plastic yeilding
                % Compute second invariant for a purely elastic stress build-up
                sxxnewe=MSXX(ym,xm)+2*MMU(MI(ym,xm))*timestep*MEXX(ym,xm)*MRAT(ym,xm);
                sxynewe=MSXY(ym,xm)+2*MMU(MI(ym,xm))*timestep*MEXY(ym,xm)*MRAT(ym,xm);
                siinewe=(sxxnewe^2+sxynewe^2)^0.5;
                % Checking yeilding criterion for strain weakening/hardening
                % C=C0, FI=FI0 for strain<=GAM0
                % C=C0+(C1-C0)/(GAM1-GAM0)*(strain-GAM0), FI=FI0+(FI1-FI0)/(GAM1-GAM0)*(strain-GAM0) for GAM0<strain<GAM1
                % C=C1, FI=FI1 for strain>=GAM0
                MCOHES=MPL(MI(ym,xm),1);
                MFRICT=MPL(MI(ym,xm),3);
                if (MGII(ym,xm)>=MPL(MI(ym,xm),6))
                    MCOHES=MPL(MI(ym,xm),2);
                    MFRICT=MPL(MI(ym,xm),4);
                end
                if (MGII(ym,xm)>MPL(MI(ym,xm),5) && MGII(ym,xm)<MPL(MI(ym,xm),6))
                    MCOHES=MPL(MI(ym,xm),1)+(MPL(MI(ym,xm),2)-MPL(MI(ym,xm),1))/(MPL(MI(ym,xm),6)-MPL(MI(ym,xm),5))*(MGII(ym,xm)-MPL(MI(ym,xm),5));
                    MFRICT=MPL(MI(ym,xm),3)+(MPL(MI(ym,xm),4)-MPL(MI(ym,xm),3))/(MPL(MI(ym,xm),6)-MPL(MI(ym,xm),5))*(MGII(ym,xm)-MPL(MI(ym,xm),5));
                end
                % Boundary  friction for sand
                if(MX(ym,xm)<0.002 || MX(ym,xm)>xsize-0.012+bintern(4)*timesum || MY(ym,xm)>ysize-0.002)
                    MCOHES=0;
                    MFRICT=sin(19/180*pi);
                end
                % Computing yelding stress for the marker
                siiyeld=MCOHES+MFRICT*MPR(ym,xm);
                if (siiyeld<0) 
                    siiyeld=0;
                end
                % Correcting rock properties for yeilding 
                if (siiyeld<siinewe && ntimestep>1 && MI(ym,xm)~=1 && MI(ym,xm)~=5)
                    % Bringing marker viscosity to yeilding stress
                    METACURP=MMU(MI(ym,xm))*timestep*siiyeld/(siinewe-siiyeld);
                    METACURP=METACURP^(1-etawt)*META(ym,xm)^etawt;
                    if(METACURP<METACUR)
                        METACUR=METACURP;
                        % Limiting viscosity for the yeilding
                        if (METACUR<etamin) 
                            METACUR=etamin;
                        end
                        if (METACUR>etamax) 
                            METACUR=etamax;
                        end
                        % Mark that plastic yeildind occur
                        plastyn=1;
                    end
                end

                % Save marker viscosity
                META(ym,xm)=METACUR;

                % Compute 1/MU values (MU is shear modulus) 
                MMUCUR=1/MMU(MI(ym,xm));


                % Add viscosity etas(), shear stress sxy(),shear modulus mus() and rock type typ() to 4 surrounding basic nodes
                % only using markers located at <=0.5 gridstep distances from nodes
                if(dx<=0.5 && dy<=0.5)
                    etas1(yn,xn)=etas1(yn,xn)+(1.0-dx).*(1.0-dy).*METACUR*mwt;
                    mus1(yn,xn)=mus1(yn,xn)+(1.0-dx).*(1.0-dy).*MMUCUR*mwt;
                    sxy1(yn,xn)=sxy1(yn,xn)+(1.0-dx).*(1.0-dy).*MSXY(ym,xm)*mwt;
                    wtetas(yn,xn)=wtetas(yn,xn)+(1.0-dx).*(1.0-dy)*mwt;
                end
                if(dx<=0.5 && dy>=0.5)
                    etas1(yn+1,xn)=etas1(yn+1,xn)+(1.0-dx).*dy.*METACUR*mwt;
                    mus1(yn+1,xn)=mus1(yn+1,xn)+(1.0-dx).*dy.*MMUCUR*mwt;
                    sxy1(yn+1,xn)=sxy1(yn+1,xn)+(1.0-dx).*dy.*MSXY(ym,xm)*mwt;
                    wtetas(yn+1,xn)=wtetas(yn+1,xn)+(1.0-dx).*dy*mwt;
                end
                if(dx>=0.5 && dy<=0.5)
                    etas1(yn,xn+1)=etas1(yn,xn+1)+dx.*(1.0-dy).*METACUR*mwt;
                    mus1(yn,xn+1)=mus1(yn,xn+1)+dx.*(1.0-dy).*MMUCUR*mwt;
                    sxy1(yn,xn+1)=sxy1(yn,xn+1)+dx.*(1.0-dy).*MSXY(ym,xm)*mwt;
                    wtetas(yn,xn+1)=wtetas(yn,xn+1)+dx.*(1.0-dy)*mwt;
                end
                if(dx>=0.5 && dy>=0.5)
                    etas1(yn+1,xn+1)=etas1(yn+1,xn+1)+dx.*dy.*METACUR*mwt;
                    mus1(yn+1,xn+1)=mus1(yn+1,xn+1)+dx.*dy.*MMUCUR*mwt;
                    sxy1(yn+1,xn+1)=sxy1(yn+1,xn+1)+dx.*dy.*MSXY(ym,xm)*mwt;
                    wtetas(yn+1,xn+1)=wtetas(yn+1,xn+1)+dx.*dy*mwt;
                end

                % Add viscosity etan(), normal stress sxx() and shear modulus mun() to the center of current cell
                etan1(yn,xn)=etan1(yn,xn)+(1.0-abs(0.5-dx)).*(1.0-abs(0.5-dy)).*METACUR*mwt;
                mun1(yn,xn)=mun1(yn,xn)+(1.0-abs(0.5-dx)).*(1.0-abs(0.5-dy)).*MMUCUR*mwt;
                sxx1(yn,xn)=sxx1(yn,xn)+(1.0-abs(0.5-dx)).*(1.0-abs(0.5-dy)).*MSXX(ym,xm)*mwt;
                wtetan(yn,xn)=wtetan(yn,xn)+(1.0-abs(0.5-dx)).*(1.0-abs(0.5-dy))*mwt;
            end

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
                mus1(i,j)=1/(mus1(i,j)./wtetas(i,j));
                sxy1(i,j)=sxy1(i,j)./wtetas(i,j);
            else
                % If no new value is interpolated from markers old value is used
                etas1(i,j)=etas0(i,j);
                mus1(i,j)=mus0(i,j);
                sxy1(i,j)=sxy0(i,j);
            end
            % Normal viscosity
            if (i<ynum && j<xnum)
                if (wtetan(i,j)~=0)
                    % Compute new value interpolated from markers
                    etan1(i,j)=etan1(i,j)./wtetan(i,j);
                    mun1(i,j)=1/(mun1(i,j)./wtetan(i,j));
                    sxx1(i,j)=sxx1(i,j)./wtetan(i,j);
                else
                    % If no new value is interpolated from markers old value is used
                    etan1(i,j)=etan0(i,j);
                    mun1(i,j)=mun0(i,j);
                    sxx1(i,j)=sxx0(i,j);
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

    % Computing viscoelastic (numerical) viscosity and stress
    % Shear stress
    for i=1:1:ynum;
        for j=1:1:xnum;
            %Viscoelasticity factor
            xelvis=etas1(i,j)/(etas1(i,j)+timestep*mus1(i,j));
            % Viscoelastic viscosity = (1-xelvis)*ETA
            etas0(i,j)=etas1(i,j)*(1-xelvis);
            % Vsicoelastic stress = xelvis*Sxy
            sxy0(i,j)=sxy1(i,j)*xelvis;
        end
    end
    % Normal stress
    for i=1:1:ynum-1;
        for j=1:1:xnum-1;
            %Viscoelasticity factor
            xelvis=etan1(i,j)/(etan1(i,j)+timestep*mun1(i,j));
            % Viscoelastic viscosity = (1-xelvis)*ETA
            etan0(i,j)=etan1(i,j)*(1-xelvis);
            % Vsicoelastic stress = xelvis*Sxx
            sxx0(i,j)=sxx1(i,j)*xelvis;
        end
    end

    % Computing right part of mechanical viscoelastic equations
    % x-Stokes
    RX1=zeros(ynum+1,xnum);
    % y-Stokes
    RY1=zeros(ynum,xnum+1);
    % continuity
    RC1=zeros(ynum-1,xnum-1);
    % Grid points cycle
    for i=2:1:ynum;
        for j=2:1:xnum;
            % Right part of x-Stokes Equation
            if(j<xnum)
                RX1(i,j)=-gx*(rho1(i,j)+rho1(i-1,j))/2;
                % Adding xelvis*dSxx0/dx
                RX1(i,j)=RX1(i,j)-(sxx0(i-1,j)-sxx0(i-1,j-1))/xstpc1(j);
                % Adding xelvis*dSxy0/dy
                RX1(i,j)=RX1(i,j)-(sxy0(i,j)-sxy0(i-1,j))/ystp1(i-1);
            end
            % Right part of y-Stokes Equation
            if(i<ynum)
                RY1(i,j)=-gy*(rho1(i,j)+rho1(i,j-1))/2;
                % Adding xelvis*dSyy0/dy using that Syy0=-Sxx0 (deviatoric stress)
                RY1(i,j)=RY1(i,j)+(sxx0(i,j-1)-sxx0(i-1,j-1))/ystpc1(i);
                % Adding xelvis*dSyx0/dx using that Syx0=Sxy0 
                RY1(i,j)=RY1(i,j)-(sxy0(i,j)-sxy0(i,j-1))/xstp1(j-1);
            end
        end
    end


    % Computing velocity field
    if (movemod==0)
        % Solving of Stokes and Continuity equations on nodes
        % and computing residuals
        % by calling function Stokes_Continuity_solver_grid() 
        % with viscoelastic numerical viscosity
        % and modified right parts
        [L,R,vx1,resx1,vy1,resy1,pr1,resc1]=Stokes_Continuity_solver_sandbox(L,R,prfirst,etas0,etan0,xnum,ynum,gridx,gridy,RX1,RY1,RC1,bleft,bright,btop,bbottom,bintern);
    end
    % Solid body rotation
    if (movemod==1)
        for i=1:1:ynum+1;
            for j=1:1:xnum+1;
                % Vx
                if(j<xnum+1)
                    % Relative distance of vx node from the model center
                    dx=((j-1)*xstp-xsize/2)/(xsize/2);
                    dy=((i-1.5)*ystp-ysize/2)/(xsize/2);
                    dr=(dx^2+dy^2)^0.5;
                    % Set vx
                    vx1(i,j)=-vyright*dy;
                end
                % Vy
                if(i<ynum+1)
                    % Relative distance of vy node from the model center
                    dx=((j-1.5)*xstp-xsize/2)/(xsize/2);
                    dy=((i-1)*ystp-ysize/2)/(xsize/2);
                    dr=(dx^2+dy^2)^0.5;
                    % Set vy
                    vy1(i,j)=vyright*dx;
                end
            end
        end
    end

    % Computing EPS'xx=-EPS'yy, EPSxy=EPSyx deviatoric strain rate tensor components from vx, vy
    % Computing spin tensor Espin
    exy = zeros(ynum,xnum);
    exx = zeros(ynum-1,xnum-1);
    esp = zeros(ynum,xnum);
    eii = zeros(ynum-1,xnum-1);
    % Grid points cycle
    for i=1:1:ynum;
        for j=1:1:xnum;
            % EPS'xx=-EPS'yy=1/2(dvx/dx-dvy/dy)
            if(i<ynum && j<xnum)
                exx(i,j)=0.5*((vx1(i+1,j+1)-vx1(i+1,j))/xstp1(j)-(vy1(i+1,j+1)-vy1(i,j+1))/ystp1(i));
            end
            % EPSxy=EPSyx=1/2(dvx/dy+dvy/dx)
            exy(i,j)=0.5*((vx1(i+1,j)-vx1(i,j))/ystpc1(i)+(vy1(i,j+1)-vy1(i,j))/xstpc1(j));
            % Espin=1/2(dvy/dx-dvx/dy) i.e. positive for clockwise rotation
            % (when x axis is directed rightward and y axis is directed downward) 
            esp(i,j)=0.5*((vy1(i,j+1)-vy1(i,j))/xstpc1(j)-(vx1(i+1,j)-vx1(i,j))/ystpc1(i));
            % EPSii=(EPS'xx^2+EPSxy^2)^0.5
            if(i>1 && j>1)
                eii(i-1,j-1)=(exx(i-1,j-1)^2+(exy(i-1,j-1)^2+exy(i,j-1)^2+exy(i-1,j)^2+exy(i,j)^2)/4)^0.5;
            end
        end
    end

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
    % Defining displacement timestep
    timestep=timestep % final displacement step
    timestepd=timestep;
    
    % Computing new stresses and stress change using the displacement timestep
    sxy2 = zeros(ynum,xnum);
    sxx2 = zeros(ynum-1,xnum-1);
    % Shear stress
    for i=1:1:ynum;
        for j=1:1:xnum;
            %Viscoelasticity factor
            xelvis=etas1(i,j)/(etas1(i,j)+timestep*mus1(i,j));
            % New viscoelastic stress = (1-xelvis)*2*ETA*EPSxy + xelvis*Sxy0
            sxy2(i,j)=(1-xelvis)*2*etas1(i,j)*exy(i,j)+xelvis*sxy1(i,j);
        end
    end
    % Normal stress
    for i=1:1:ynum-1;
        for j=1:1:xnum-1;
            %Viscoelasticity factor
            xelvis=etan1(i,j)/(etan1(i,j)+timestep*mun1(i,j));
            % New viscoelastic stress = (1-xelvis)*2*ETA*EPSxx + xelvis*Sxx0
            sxx2(i,j)=(1-xelvis)*2*etan1(i,j)*exx(i,j)+xelvis*sxx1(i,j);
        end
    end
    % Stress change
    dsxy = sxy2-sxy1;
    dsxx = sxx2-sxx1;
    
    % Computing strain rate and pressure for markers
    for xm = 1:1:mxnum
        for ym = 1:1:mynum  

            % Check markers inside the grid
            if (MX(ym,xm)>=0 && MX(ym,xm)<=xsize && MY(ym,xm)>=0 && MY(ym,xm)<=ysize) 

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

                % Define indexes for upper left BASIC node in the cell where the marker is
                xnmin=MXN(ym,xm);
                ynmin=MYN(ym,xm);

                % Calculating strain rate for marker
                %
                % Interpolating squares of EPS'xx=-EPS'yy from cell centers
                % EPS'xx-nodes are displaced rightward and downward for 1/2 of gridsteps
                % Horizontal EPS'xx index
                xn=xnmin;
                if(MX(ym,xm)<gridcx(xn+1))
                    xn=xn-1;
                end
                if (xn<1)
                    xn=1;
                end
                if (xn>xnum-2)
                    xn=xnum-2;
                end
                % Vertical EPS'xx index
                yn=ynmin;
                if(MY(ym,xm)<gridcy(yn+1))
                    yn=yn-1;
                end
                if (yn<1)
                    yn=1;
                end
                if (yn>ynum-2)
                    yn=ynum-2;
                end

                % Define and check normalized distances from marker to the upper left EPS'xx-node;
                dx=(MX(ym,xm)-gridcx(xn+1))./xstpc1(xn+1);
                dy=(MY(ym,xm)-gridcy(yn+1))./ystpc1(yn+1);

                % Calculate and save Marker EPS'xx from four surrounding nodes
                exxm=0;
                exxm=exxm+(1.0-dx).*(1.0-dy).*exx(yn,xn);
                exxm=exxm+(1.0-dx).*dy.*exx(yn+1,xn);
                exxm=exxm+dx.*(1.0-dy).*exx(yn,xn+1);
                exxm=exxm+dx.*dy.*exx(yn+1,xn+1);
                MEXX(ym,xm)=exxm;
                % Calculate Marker SIG'xx from four surrounding nodes
                sxxm=0;
                sxxm=sxxm+(1.0-dx).*(1.0-dy).*sxx2(yn,xn);
                sxxm=sxxm+(1.0-dx).*dy.*sxx2(yn+1,xn);
                sxxm=sxxm+dx.*(1.0-dy).*sxx2(yn,xn+1);
                sxxm=sxxm+dx.*dy.*sxx2(yn+1,xn+1);
                % Calculate Marker dSIG'xx from four surrounding nodes
                dsxxm=0;
                dsxxm=dsxxm+(1.0-dx).*(1.0-dy).*dsxx(yn,xn);
                dsxxm=dsxxm+(1.0-dx).*dy.*dsxx(yn+1,xn);
                dsxxm=dsxxm+dx.*(1.0-dy).*dsxx(yn,xn+1);
                dsxxm=dsxxm+dx.*dy.*dsxx(yn+1,xn+1);
                

                % Calculate and save Marker pressure from four surrounding nodes
                prm=0;
                prm=prm+(1.0-dx).*(1.0-dy).*pr1(yn,xn);
                prm=prm+(1.0-dx).*dy.*pr1(yn+1,xn);
                prm=prm+dx.*(1.0-dy).*pr1(yn,xn+1);
                prm=prm+dx.*dy.*pr1(yn+1,xn+1);
                MPR(ym,xm)=prm;


                % Interpolating EPSxy=EPSyx from basic nodes
                % Horizontal EPSxy index
                xn=xnmin;
                % Vertical EPSxy index
                yn=ynmin;

                % Define and check normalized distances from marker to the upper left VX-node;
                dx=(MX(ym,xm)-gridx(xn))./xstp1(xn);
                dy=(MY(ym,xm)-gridy(yn))./ystp1(yn);

                % Calculate and save Marker EPSxy from four surrounding nodes
                exym=0;
                exym=exym+(1.0-dx).*(1.0-dy).*exy(yn,xn);
                exym=exym+(1.0-dx).*dy.*exy(yn+1,xn);
                exym=exym+dx.*(1.0-dy).*exy(yn,xn+1);
                exym=exym+dx.*dy.*exy(yn+1,xn+1);
                MEXY(ym,xm)=exym;
                % Calculate Marker SIGxy from four surrounding nodes
                sxym=0;
                sxym=sxym+(1.0-dx).*(1.0-dy).*sxy2(yn,xn);
                sxym=sxym+(1.0-dx).*dy.*sxy2(yn+1,xn);
                sxym=sxym+dx.*(1.0-dy).*sxy2(yn,xn+1);
                sxym=sxym+dx.*dy.*sxy2(yn+1,xn+1);
                % Calculate Marker SIGxy from four surrounding nodes
                dsxym=0;
                dsxym=dsxym+(1.0-dx).*(1.0-dy).*dsxy(yn,xn);
                dsxym=dsxym+(1.0-dx).*dy.*dsxy(yn+1,xn);
                dsxym=dsxym+dx.*(1.0-dy).*dsxy(yn,xn+1);
                dsxym=dsxym+dx.*dy.*dsxy(yn+1,xn+1);
                
                % Computing second strain rate invariant 
                % for the marker using grid values
                eiimg=(MEXX(ym,xm)^2+MEXY(ym,xm)^2)^0.5;
                % Correcting strain rate for the marker using Maxwell model  
                if (eiimg>0)
                    % Computing second strain rate invariant for the marker
                    % from stresses using Maxwell model
                    eiim=((sxxm/2/META(ym,xm)+dsxxm/2/timestep/MMU(MI(ym,xm)))^2+(sxym/2/META(ym,xm)+dsxym/2/timestep/MMU(MI(ym,xm)))^2)^0.5;
                    % Computing EiiMarker/EiiGrid ratio
                    MRAT(ym,xm)=(eiim/eiimg);
                else
                    MRAT(ym,xm)=1;
                end
            end
        end
    end


    % Computing subgrid stress changes for markers
    if (dsubgrids>0)
        % Clear subgrid stress changes for nodes
        dsxyn=zeros(ynum,xnum);
        dsxxn=zeros(ynum-1,xnum-1);
        % Clear wights for Sxy
        wtetas=zeros(ynum,xnum);
        % Clear wights for Sxx
        wtetan=zeros(ynum-1,xnum-1);
        % Marker cycle
        for xm = 1:1:mxnum
            for ym = 1:1:mynum  

            % Check markers inside the grid
                if (MX(ym,xm)>=0 && MX(ym,xm)<=xsize && MY(ym,xm)>=0 && MY(ym,xm)<=ysize)
                    
                    % Compute local stress relaxation timescale (Maxwell time) for the marker
                    sdm=META(ym,xm)/MMU(MI(ym,xm));
                    % Computing degree of subgrid stress relaxation
                    sdif=-dsubgrids*timestep/sdm;
                    if(sdif<-30) 
                        sdif=-30;
                    end
                    sdif=(1-exp(sdif));

                    %  yn    sxy(yn,xn)--------------------sxy(yn,xn+1)
                    %           ?           ^                  ?
                    %           ?           ?                  ?
                    %           ?          dy                  ?
                    %           ?           ?                  ?
                    %           ?           v                  ?
                    %           ?<----dx--->o MSXY(ym,xm)      ?
                    %           ?                              ?
                    %           ?                              ?
                    %  yn+1  sxy(yn+1,xn)-------------------sxy(yn+1,xn+1)
                    %
                    %
                    % Interpolating old shear stress from Sxy nodes
                    %
                    % Define indexes for upper left node in the cell where the marker is
                    xn=MXN(ym,xm);
                    yn=MYN(ym,xm);

                    % Define normalized distances from marker to the upper left node;
                    dx=(MX(ym,xm)-gridx(xn))./xstp1(xn);
                    dy=(MY(ym,xm)-gridy(yn))./ystp1(yn);

                    % Compute marker weight koefficient from cell dimensions
                    % Number of markers in a cell is in invert proportion to the cell volume
                    mwt=1;%/xstp1(xn)/ystp1(yn);


                    % Interpolate old Sxy stress for the marker
                    sxym=0;
                    sxym=sxym+(1.0-dx).*(1.0-dy).*sxy1(yn,xn);
                    sxym=sxym+(1.0-dx).*dy.*sxy1(yn+1,xn);
                    sxym=sxym+dx.*(1.0-dy).*sxy1(yn,xn+1);
                    sxym=sxym+dx.*dy.*sxy1(yn+1,xn+1);
                    % Calculate Nodal-Marker subgrid Sxy stress difference
                    dsxym=sxym-MSXY(ym,xm);
                    % Relaxing Nodal-Marker subgrid Sxy stress difference
                    dsxym=dsxym*sdif;    

                    % Correcting old stress for the marker
                    MSXY(ym,xm)=MSXY(ym,xm)+dsxym;

                    % Interpolating subgrid Sxy stress changes to 4 nodes
                    % only using markers located at <=0.5 gridstep distances from nodes
                    if(dx<=0.5 && dy<=0.5)
                        dsxyn(yn,xn)=dsxyn(yn,xn)+(1.0-dx).*(1.0-dy).*dsxym*mwt;
                        wtetas(yn,xn)=wtetas(yn,xn)+(1.0-dx).*(1.0-dy)*mwt;
                    end
                    if(dx<=0.5 && dy>=0.5)
                        dsxyn(yn+1,xn)=dsxyn(yn+1,xn)+(1.0-dx).*dy.*dsxym*mwt;
                        wtetas(yn+1,xn)=wtetas(yn+1,xn)+(1.0-dx).*dy*mwt;
                    end
                    if(dx>=0.5 && dy<=0.5)
                        dsxyn(yn,xn+1)=dsxyn(yn,xn+1)+dx.*(1.0-dy).*dsxym*mwt;
                        wtetas(yn,xn+1)=wtetas(yn,xn+1)+dx.*(1.0-dy)*mwt;
                    end
                    if(dx>=0.5 && dy>=0.5)
                        dsxyn(yn+1,xn+1)=dsxyn(yn+1,xn+1)+dx.*dy.*dsxym*mwt;
                        wtetas(yn+1,xn+1)=wtetas(yn+1,xn+1)+dx.*dy*mwt;
                    end
                    
                    % Computing marker weight for the center of current
                    % basic cell where Sxx stress is located
                    mwt=mwt*(1.0-abs(0.5-dx)).*(1.0-abs(0.5-dy));
                    %  yn    sxx(yn,xn)--------------------sxx(yn,xn+1)
                    %           ?           ^                  ?
                    %           ?           ?                  ?
                    %           ?          dy                  ?
                    %           ?           ?                  ?
                    %           ?           v                  ?
                    %           ?<----dx--->o MSXX(ym,xm)       ?
                    %           ?                              ?
                    %           ?                              ?
                    %  yn+1  sxx(yn+1,xn)-------------------sxx(yn+1,xn+1)
                    %
                    %
                    % Interpolating old normal stress from Sxx nodes
                    %
                    % Define, check indexes for upper left node in the Sxx cell where the marker is
                    if (MX(ym,xm)<gridcx(xn+1))
                        xn=xn-1;
                    end
                    if(xn<1)
                        xn=1;
                    end
                    if(xn>xnum-2)
                        xn=xnum-2;
                    end
                    if (MY(ym,xm)<gridcy(yn+1))
                        yn=yn-1;
                    end
                    if(yn<1)
                        yn=1;
                    end
                    if(yn>ynum-2)
                        yn=ynum-2;
                    end

                    % Define normalized distances from marker to the upper left node;
                    dx=(MX(ym,xm)-gridcx(xn+1))./xstpc1(xn+1);
                    dy=(MY(ym,xm)-gridcy(yn+1))./ystpc1(yn+1);

                    % Interpolate old Sxx stress for the marker
                    sxxm=0;
                    sxxm=sxxm+(1.0-dx).*(1.0-dy).*sxx1(yn,xn);
                    sxxm=sxxm+(1.0-dx).*dy.*sxx1(yn+1,xn);
                    sxxm=sxxm+dx.*(1.0-dy).*sxx1(yn,xn+1);
                    sxxm=sxxm+dx.*dy.*sxx1(yn+1,xn+1);
                    % Calculate Nodal-Marker subgrid Sxx stress difference
                    dsxxm=sxxm-MSXX(ym,xm);
                    % Relaxing Nodal-Marker subgrid Sxx stress difference
                    dsxxm=dsxxm*sdif;    

                    % Correcting old stress for the marker
                    MSXX(ym,xm)=MSXX(ym,xm)+dsxxm;

                    % Interpolating subgrid Sxx stress changes for the center of current basic cell
                    xn=MXN(ym,xm);
                    yn=MYN(ym,xm);
                    dsxxn(yn,xn)=dsxxn(yn,xn)+dsxxm*mwt;
                    wtetan(yn,xn)=wtetan(yn,xn)+mwt;

                end
            end
        end

        % Computing subgrid stress changes for nodes
        for i=1:1:ynum;
            for j=1:1:xnum;
                % Density
                if (wtetas(i,j)~=0)
                    % Compute new value interpolated from markers
                    dsxyn(i,j)=dsxyn(i,j)./wtetas(i,j);
                end
                if (j<xnum && i<ynum && wtetan(i,j)~=0)
                    % Compute new value interpolated from markers
                    dsxxn(i,j)=dsxxn(i,j)./wtetan(i,j);
                end
            end
        end

        % Subtracting subgrid stress change part from nodal stress changes
        dsxy=dsxy-dsxyn;
        dsxx=dsxx-dsxxn;

    end

    % Updating stress for markers
    for xm = 1:1:mxnum
        for ym = 1:1:mynum  

        % Check markers inside the grid
            if (MX(ym,xm)>=0 && MX(ym,xm)<=xsize && MY(ym,xm)>=0 && MY(ym,xm)<=ysize)

                %  yn    sxy(yn,xn)--------------------sxy(yn,xn+1)
                %           ?           ^                  ?
                %           ?           ?                  ?
                %           ?          dy                  ?
                %           ?           ?                  ?
                %           ?           v                  ?
                %           ?<----dx--->o MSXY(ym,xm)       ?
                %           ?                              ?
                %           ?                              ?
                %  yn+1  sxy(yn+1,xn)-------------------sxy(yn+1,xn+1)
                %
                %
                % Interpolating old shear stress changes from Sxy nodes
                %
                % Define indexes for upper left node in the cell where the marker is
                xn=MXN(ym,xm);
                yn=MYN(ym,xm);

                % Define normalized distances from marker to the upper left node;
                dx=(MX(ym,xm)-gridx(xn))./xstp1(xn);
                dy=(MY(ym,xm)-gridy(yn))./ystp1(yn);

                % Interpolate old Sxy stress change for the marker
                dsxym=0;
                dsxym=dsxym+(1.0-dx).*(1.0-dy).*dsxy(yn,xn);
                dsxym=dsxym+(1.0-dx).*dy.*dsxy(yn+1,xn);
                dsxym=dsxym+dx.*(1.0-dy).*dsxy(yn,xn+1);
                dsxym=dsxym+dx.*dy.*dsxy(yn+1,xn+1);

                % Update stress for the marker
                MSXY(ym,xm)=MSXY(ym,xm)+dsxym;

                %  yn    sxx(yn,xn)--------------------sxx(yn,xn+1)
                %           ?           ^                  ?
                %           ?           ?                  ?
                %           ?          dy                  ?
                %           ?           ?                  ?
                %           ?           v                  ?
                %           ?<----dx--->o MSXX(ym,xm)       ?
                %           ?                              ?
                %           ?                              ?
                %  yn+1  sxx(yn+1,xn)-------------------sxx(yn+1,xn+1)
                %
                %
                % Interpolating old normal stress changes from Sxx nodes
                %
                % Define, check indexes for upper left node in the Sxx cell where the marker is
                if (MX(ym,xm)<gridcx(xn+1))
                    xn=xn-1;
                end
                if(xn<1)
                    xn=1;
                end
                if(xn>xnum-2)
                    xn=xnum-2;
                end
                if (MY(ym,xm)<gridcy(yn+1))
                    yn=yn-1;
                end
                if(yn<1)
                    yn=1;
                end
                if(yn>ynum-2)
                    yn=ynum-2;
                end

                % Define normalized distances from marker to the upper left node;
                dx=(MX(ym,xm)-gridcx(xn+1))./xstpc1(xn+1);
                dy=(MY(ym,xm)-gridcy(yn+1))./ystpc1(yn+1);

                % Interpolate old Sxx stress for the marker
                dsxxm=0;
                dsxxm=dsxxm+(1.0-dx).*(1.0-dy).*dsxx(yn,xn);
                dsxxm=dsxxm+(1.0-dx).*dy.*dsxx(yn+1,xn);
                dsxxm=dsxxm+dx.*(1.0-dy).*dsxx(yn,xn+1);
                dsxxm=dsxxm+dx.*dy.*dsxx(yn+1,xn+1);

                % Correcting old stress for the marker
                MSXX(ym,xm)=MSXX(ym,xm)+dsxxm;

            end
        end
    end


    
    % Solving Temperature equation
    if (timestep>0 && tempmax>0)
        
        % Computing right part of temperature equation
        RT1=zeros(ynum,xnum);
        % Computing viscoelastic shear heating for Temperature nodes
        % Hs=2*Sxx*Sxx/2/etan+2*Sxy*Sxy/2/etas
        % Grid points cycle
        for i=2:1:ynum-1;
            for j=2:1:xnum-1;
                % Shear heating on(1)/off(0)
                if(frictyn==1)
                    % Adding 2*Sxy*Sxy/2/etas
                    RT1(i,j)=RT1(i,j)+sxy2(i,j)^2/etas1(i,j);
                    % Computing and adding 2*Sxx*Sxx/2/etan
                    RT1(i,j)=RT1(i,j)+(sxx2(i-1,j-1)^2/etan1(i-1,j-1)+sxx2(i,j-1)^2/etan1(i,j-1)+sxx2(i-1,j)^2/etan1(i-1,j)+sxx2(i,j)^2/etan1(i,j))/4;
                end
            end
        end

        % Solving temperature equation making (if needed) several thermal
        % timesteps for one displacement timestep
        % Set current thermal timestep 
        timestept=timestep
        % Set total thermal timestep
        timesteps=0;
        % Set old Temperature
        tk0=tk1;
        while (timesteps<timestep)
            % Solving Temperature equation with thermal timestep
            [LT,RT,tk2,rest]=Temperature_solver_sandbox(LT,RT,timestept,xnum,ynum,gridx,gridy,kt1,rhocp1,tk0,RT1,bleftt,brightt,btopt,bbottomt);
            % Computing temperature changes
            dtk1=tk2-tk0;
            % Checking temperature changes
            dtkmax=max(max(abs(dtk1)))
            % Repeating temperature solution if temperature changes are too big
            if(dtkmax>tempmax)
                % Computing reduced timestep
                timestept=timestept*tempmax/dtkmax
                % Solving Temperature equation with reduced timestep
                [LT,RT,tk2,rest]=Temperature_solver_sandbox(LT,RT,timestept,xnum,ynum,gridx,gridy,kt1,rhocp1,tk0,RT1,bleftt,brightt,btopt,bbottomt);
                % Computing temperature changes
            end
            % Add total thermal timestep
            timesteps=timesteps+timestept
            % Compute current thermal timestep
            if (timestept>timestep-timesteps)
                timestept=timestep-timesteps
            else
                timestept=timestept
            end
            % Update old temperature
            tk0=tk2;
        end
        % Compute temperature changes
        dtk1=tk2-tk1;


        % Computing subgrid diffusion for markers
        if (dsubgridt>0)
            % Clear subgrid temperature changes for nodes
            dtkn=zeros(ynum,xnum);
            % Clear wights for basic nodes
            wtnodes=zeros(ynum,xnum);
            % Marker cycle
            for xm = 1:1:mxnum
                for ym = 1:1:mynum  
                    
                % Check markers inside the grid
                    if (MX(ym,xm)>=0 && MX(ym,xm)<=xsize && MY(ym,xm)>=0 && MY(ym,xm)<=ysize) 
                        
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
                        %
                        % Define indexes for upper left node in the cell where the marker is
                        xn=MXN(ym,xm);
                        yn=MYN(ym,xm);

                        % Define normalized distances from marker to the upper left node;
                        dx=(MX(ym,xm)-gridx(xn))./xstp1(xn);
                        dy=(MY(ym,xm)-gridy(yn))./ystp1(yn);

                        % Compute marker weight koefficient from cell dimensions
                        % Number of markers in a cell is in invert proportion to the cell volume
                        mwt=1;%/xstp1(xn)/ystp1(yn);


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

                        % Interpolating subgrid temperature changes to 4 nodes
                        dtkn(yn,xn)=dtkn(yn,xn)+(1.0-dx).*(1.0-dy).*dtkm*mwt;
                        wtnodes(yn,xn)=wtnodes(yn,xn)+(1.0-dx).*(1.0-dy)*mwt;
                        
                        dtkn(yn+1,xn)=dtkn(yn+1,xn)+(1.0-dx).*dy.*dtkm*mwt;
                        wtnodes(yn+1,xn)=wtnodes(yn+1,xn)+(1.0-dx).*dy*mwt;
                        
                        dtkn(yn,xn+1)=dtkn(yn,xn+1)+dx.*(1.0-dy).*dtkm*mwt;
                        wtnodes(yn,xn+1)=wtnodes(yn,xn+1)+dx.*(1.0-dy)*mwt;
                        
                        dtkn(yn+1,xn+1)=dtkn(yn+1,xn+1)+dx.*dy.*dtkm*mwt;
                        wtnodes(yn+1,xn+1)=wtnodes(yn+1,xn+1)+dx.*dy*mwt;
                        
                    end
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

                % Check markers inside the grid
                if (MX(ym,xm)>=0 && MX(ym,xm)<=xsize && MY(ym,xm)>=0 && MY(ym,xm)<=ysize) 

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
                    %
                    % Define indexes for upper left node in the cell where the marker is
                    xn=MXN(ym,xm);
                    yn=MYN(ym,xm);

                    % Define normalized distances from marker to the upper left node;
                    dx=(MX(ym,xm)-gridx(xn))./xstp1(xn);
                    dy=(MY(ym,xm)-gridy(yn))./ystp1(yn);

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

    end
    
   

    % Moving Markers by velocity field
    if(markmove>0)
        % Create arrays for velocity and spin of markers
        vxm=zeros(4,1);
        vym=zeros(4,1);
        espm=zeros(4,1);
        % Marker cycle
        for xm = 1:1:mxnum
            for ym = 1:1:mynum  

                % Check markers inside the grid
                if (MX(ym,xm)>=0 && MX(ym,xm)<=xsize && MY(ym,xm)>=0 && MY(ym,xm)<=ysize) 

                    % Save marker coordinates
                    xcur=MX(ym,xm);
                    ycur=MY(ym,xm);
                    % Defining number of Runge-Kutta cycles
                    for rk=1:1:markmove

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

                       % Define indexes for upper left BASIC node in the cell where the marker is
                        % using bisection
                        % Load horizontal and vertical indexes
                        if(rk==1)
                            xnmin=MXN(ym,xm);
                            ynmin=MYN(ym,xm);
                        else
                            % Find horizontal index
                            xnmin=1;
                            xnmax=xnum;
                            while ((xnmax-xnmin)>1)
                                % !!! SUBTRACT 0.5 since int16(0.5)=1
                                xn=double(int16((xnmax+xnmin)./2-0.5));
                                if(gridx(xn)>xcur)
                                    xnmax=xn;
                                else
                                    xnmin=xn;
                                end
                            end
                            % Check horizontal index
                            if (xnmin<1)
                                xnmin=1;
                            end
                            if (xnmin>xnum-1)
                                xnmin=xnum-1;
                            end
                            % Find vertical index
                            ynmin=1;
                            ynmax=ynum;
                            while ((ynmax-ynmin)>1)
                                % !!! SUBTRACT 0.5 since int16(0.5)=1
                                yn=double(int16((ynmax+ynmin)./2-0.5));
                                if(gridy(yn)>ycur)
                                    ynmax=yn;
                                else
                                    ynmin=yn;
                                end
                            end
                            % Check vertical index
                            if (ynmin<1)
                                ynmin=1;
                            end
                            if (ynmin>ynum-1)
                                ynmin=ynum-1;
                            end
                        end

                        % Define indexes for upper left node in the Vx-cell where the marker is
                        % Horizontal Vx index
                        xn=xnmin;
                        % Vertical Vx index
                        yn=ynmin;
                        if(ycur>gridcy(yn+1))
                            yn=yn+1;
                        end
                        if (yn>ynum)
                            yn=ynum;
                        end
                        
                        % Define and check normalized distances from marker to the upper left VX-node;
                        dx=(xcur-gridx(xn))./xstp1(xn);
                        dy=(ycur-gridcy(yn))./ystpc1(yn);
                        
                        % Calculate Marker velocity from four surrounding Vx nodes
                        vxm(rk)=0;
                        vxm(rk)=vxm(rk)+(1.0-dx).*(1.0-dy).*vx1(yn,xn);
                        vxm(rk)=vxm(rk)+(1.0-dx).*dy.*vx1(yn+1,xn);
                        vxm(rk)=vxm(rk)+dx.*(1.0-dy).*vx1(yn,xn+1);
                        vxm(rk)=vxm(rk)+dx.*dy.*vx1(yn+1,xn+1);

                        % Define indexes for upper left node in the VY-cell where the marker is
                        % Vertical Vy index
                        yn=ynmin;
                        % Horizontal Vy index
                        xn=xnmin;
                        if(xcur>gridcx(xn+1))
                            xn=xn+1;
                        end
                        if (xn>xnum)
                            xn=xnum;
                        end
                        
                        % Define and check normalized distances from marker to the upper left VX-node;
                        dx=(xcur-gridcx(xn))./xstpc1(xn);
                        dy=(ycur-gridy(yn))./ystp1(yn);
                        
                        % Calculate Marker velocity from four surrounding nodes
                        vym(rk)=0;
                        vym(rk)=vym(rk)+(1.0-dx).*(1.0-dy).*vy1(yn,xn);
                        vym(rk)=vym(rk)+(1.0-dx).*dy.*vy1(yn+1,xn);
                        vym(rk)=vym(rk)+dx.*(1.0-dy).*vy1(yn,xn+1);
                        vym(rk)=vym(rk)+dx.*dy.*vy1(yn+1,xn+1);
                        
                        % Define indexes for upper left node in the Espin cell where the marker is
                        xn=xnmin;
                        yn=ynmin;

                        % Define normalized distances from marker to the upper left node;
                        dx=(xcur-gridx(xn))./xstp1(xn);
                        dy=(ycur-gridy(yn))./ystp1(yn);

                        % Interpolate old Sxy stress for the marker
                        espm(rk)=0;
                        espm(rk)=espm(rk)+(1.0-dx).*(1.0-dy).*esp(yn,xn);
                        espm(rk)=espm(rk)+(1.0-dx).*dy.*esp(yn+1,xn);
                        espm(rk)=espm(rk)+dx.*(1.0-dy).*esp(yn,xn+1);
                        espm(rk)=espm(rk)+dx.*dy.*esp(yn+1,xn+1);
                       

                        % Update coordinates for the next cycle
                        if(rk<4)
                            if (rk<3)
                                xcur=MX(ym,xm)+timestep/2*vxm(rk);
                                ycur=MY(ym,xm)+timestep/2*vym(rk);
                            else
                                xcur=MX(ym,xm)+timestep*vxm(rk);
                                ycur=MY(ym,xm)+timestep*vym(rk);
                            end
                        end
                                              
                        
                    end
                    % Recompute velocity and spin using 4-th order Runge_Kutta
                    if (markmove==4)
                        vxm(1)=(vxm(1)+2*vxm(2)+2*vxm(3)+vxm(4))/6;
                        vym(1)=(vym(1)+2*vym(2)+2*vym(3)+vym(4))/6;
                        espm(1)=(espm(1)+2*espm(2)+2*espm(3)+espm(4))/6;
                    end

                    % Displacing Marker according to its velocity
                    MX(ym,xm)=MX(ym,xm)+timestep*vxm(1);
                    MY(ym,xm)=MY(ym,xm)+timestep*vym(1);
                    
                    % Rotate stress on marker according to its spin
                    % Compute amount of rotation from spin rate:
                    % Espin=1/2(dvy/dx-dvx/dy) i.e. positive for clockwise rotation
                    % (when x axis is directed rightward and y axis is directed downward) 
                    espm(1)=espm(1)*timestep;
                    % Save old stresses
                    msxxold=MSXX(ym,xm);
                    msxyold=MSXY(ym,xm);
                    % SxyNEW=0.5(Sxx-Syy)*sin(2*Espin*dt)+Sxy*cos(2*Espin*dt)
                    % where Sxx-Syy=2Sxx
                    MSXY(ym,xm)=msxxold*sin(2*espm(1))+msxyold*cos(2*espm(1));
                    % SxxNEW=Sxx*(cos(Espin*dt))^2+Syy*(sin(Espin*dt))^2-Sxy*sin(2*Espin*dt)
                    % where Sxx=-Syy
                    MSXX(ym,xm)=msxxold*((cos(espm(1)))^2-(sin(espm(1)))^2)-msxyold*sin(2*espm(1));
                    
                    % Adding marker strain based on grid strain rates
                    MGII(ym,xm)=MGII(ym,xm)+timestep*(MEXX(ym,xm)^2+MEXY(ym,xm)^2)^0.5;
   
                end
                
           end
        end
    end


    % Interpolating vx, vy for the basic grid
    vxb=zeros(ynum,xnum);
    vyb=zeros(ynum,xnum);
    for j=1:1:xnum
        for i=1:1:ynum
            vxb(i,j)=(vx1(i,j)+vx1(i+1,j))/2;
            vyb(i,j)=(vy1(i,j)+vy1(i,j+1))/2;
        end
    end
    
   
    figure(1);colormap('Jet');clf
  
    % Plotting Sxx
    subplot(4,1,1)
    pcolor(sxx2);
    shading interp;
    axis tight;
    title('Sxx')
    colorbar; 
    axis ij image;

    % Plotting Sxy
    subplot(4,1,2)
    pcolor(sxy2);
    shading interp;
    axis tight;
    title('Sxy')
    colorbar;    
    axis ij image;

    % Plotting Pr
    subplot(4,1,3)
    pcolor(pr1);
    shading interp;
    axis tight;
    title('Pr')
    colorbar;    
    axis ij image;
    
    % Plotting Eii
    subplot(4,1,4)
    pcolor(log10(eii));
    shading interp;
    axis tight;
    title('Eii')
    colorbar;    
    axis ij image;
    
    
    figure(2);
    subplot(2,1,1);colormap('Jet');clf
    % Plotting viscoelastic viscosity
%    pcolor(gridcx(2:1:xnum)/1000,gridcy(2:1:ynum)/1000,sxx2);
    pcolor(gridx/1000,gridy/1000,log10(etas0));
    shading interp;
    colorbar;
    title(['Step=',num2str(ntimestep),' Myr=',num2str(timesum*1e-6/(365.25*24*3600))]); 
    hold on;
    % Plotting velocity
    quiver(gridx(2:4:xnum-1)/1000,gridy(2:4:ynum-1)/1000,vxb(2:4:ynum-1,2:4:xnum-1),vyb(2:4:ynum-1,2:4:xnum-1),'k');
    hold off;
    axis ij image;
    axis([0 xsize/1000 0 ysize/1000]);
    box on;
    subplot(2,1,2);
    % Plotting viscosity
    pcolor(gridx/1000,gridy/1000,log10(etas1));
    shading interp;
    colorbar;
    axis ij image;
    axis([0 xsize/1000 0 ysize/1000]);
    box on;
%     print ('-dtiff', '-r300','-zbuffer ','fig3');

    
    % Visualizing marker type
    % Pixel grid resolution
    xresol=mxnum+1;
    yresol=mynum+1;
    ngrid=2;
    sxstp=xsize/(xresol-1);
    systp=ysize/(yresol-1);
    mx=0:sxstp*100:xsize*100;
    my=0:systp*100:ysize*100;
    % Process markers
    markcom=NaN*ones(yresol,xresol);
    markdis=1e+20*ones(yresol,xresol);
    markgii=NaN*ones(yresol,xresol);
    for xm = 1:1:mxnum
        for ym = 1:1:mynum 
            % Define pixel cell
            m1=double(int16(MX(ym,xm)/sxstp-0.5))+1;
            m2=double(int16(MY(ym,xm)/systp-0.5))+1;
            if (m1<1)
                m1=1;
            end
            if (m1>xresol-1)
                m1=xresol-1;
            end
            if (m2<1)
                m2=1;
            end
            if (m2>yresol-1)
                m2=yresol-1;
            end
            % Define indexes of surrounding pixels
            m10min=m1-ngrid;
            if (m10min<1)
                m10min=1;
            end
            m10max=m1+1+ngrid;
            if (m10max>xresol)
                m10max=xresol;
            end
            m20min=m2-ngrid;
            if (m20min<1)
                m20min=1;
            end
            m20max=m2+1+ngrid;
            if (m20max>yresol)
                m20max=yresol;
            end
            % Update pixels around the marker
            for m10 = m10min:1:m10max
                for m20 = m20min:1:m20max 
                    % Check distance to current pixel
                    dx=MX(ym,xm)-(m10-1)*sxstp;
                    dy=MY(ym,xm)-(m20-1)*systp;
                    dd=(dx*dx+dy*dy)^0.5;
                    if(dd<markdis(m20,m10))
                        markcom(m20,m10)=MI(ym,xm);
                        markgii(m20,m10)=MGII(ym,xm);
                        markdis(m20,m10)=dd;
                    end
                end
            end
         end
    end
    
    % Draw composition
    figure(3);colormap(gray);clf
    pcolor(mx,my,-markcom);
    box on;
    shading interp;
    axis ij image;
    title(['Step=',num2str(ntimestep),' Myr=',num2str(timesum*1e-6/(365.25*24*3600))]); 

    % Draw strain
    figure(4);colormap(gray);clf
    pcolor(mx,my,-log10(markgii));
    box on;
    colorbar;
    shading interp;
    axis ij image;

    
    
    timesum=timesum
    timestep=timestep
    
    % Advance in time
    timesum=timesum+timestep
    
    ntimestep=ntimestep

       
    
    % Increase timestep
    if(ntimestep==10 || ntimestep==20 || ntimestep==30 || ntimestep==40)
        timemax=timemax*10^0.5;
        markmax=markmax*10^0.5;
    end
    pause(0.1)
    
    % Overal distance
    dxshort=-bintern(4)*timesum
    if(dxshort>0.02)
        break;
    end
end