% Solving Navier-Stokes, continuity and temperature eq.
% in primitive variable formulation
% with variable viscosity and thermal conductivity 
% for visco-elasto-plastic material
% using  FD with staggered grid
% and marker in cell techique

% Clearing memory and figures
clear all; clf

% Define Numerical model
xsize=100000; % Horizontal model size, m
ysize=100000; % Vertical model size, m
Nx=101; % Horizontal grid resolution
Ny=101; % Vertical grid resolution
Nx1=Nx+1;
Ny1=Ny+1;
dx=xsize/(Nx-1); % Horizontal grid step, m
dy=ysize/(Ny-1); % Vertical grid step, m

% Define Gravity
gx=0; % Horizontal gravity acceleration, m/s^2
gy=0; % Vertical gravity acceleration, m/s^2

% Coordinates of different nodal points
% Basic nodes
x=0:dx:xsize; % Horizontal coordinates of basic grid points, m
y=0:dy:ysize; % Vertical coordinates of basic grid points, m
% Vx-Nodes
xvx=0:dx:xsize+dy; % Horizontal coordinates of vx grid points, m
yvx=-dy/2:dy:ysize+dy/2; % Vertical coordinates of vx grid points, m
% Vy-nodes
xvy=-dx/2:dx:xsize+dx/2; % Horizontal coordinates of vy grid points, m
yvy=0:dy:ysize+dy; % Vertical coordinates of vy grid points, m
% P-Nodes
xp=-dx/2:dx:xsize+dx/2; % Horizontal coordinates of P grid points, m
yp=-dy/2:dy:ysize+dy/2; % Vertical coordinates of P grid points, m

% Nodal arrays
% Basic nodes
ETA=zeros(Ny,Nx); % Viscoplastic Viscosity, Pa*s
ETA0=zeros(Ny,Nx); % Viscous Viscosity, Pa*s
GGG=zeros(Ny,Nx); % Shear modulus, Pa
EXY=zeros(Ny,Nx); % EPSILONxy, 1/s
SXY=zeros(Ny,Nx); % SIGMAxy, 1/s
SXY0=zeros(Ny,Nx); % SIGMA0xy, 1/s
wyx=zeros(Ny,Nx); % Rotation rate, 1/s
COH=zeros(Ny,Nx); % Cohesion, Pa
FRIS=zeros(Ny,Nx); % Static Friction
FRID=zeros(Ny,Nx); % Dynamic Friction
YNY=zeros(Ny,Nx); % Plastic yielding mark, 1=yes,0=no
% Vx-Nodes
RHOX=zeros(Ny1,Nx1); % Density, kg/m^3
KX=zeros(Ny1,Nx1); % Thermal conductivity, W/m/K
vx=zeros(Ny1,Nx1); % vx-velocity m/s
vx0=zeros(Ny1,Nx1); % Old vx-velocity m/s
% Vy-Nodes
RHOY=zeros(Ny1,Nx1); % Density, kg/m^3
KY=zeros(Ny1,Nx1); % Thermal conductivity, W/m/K
vy=zeros(Ny1,Nx1); % vy-velocity m/s
vy0=zeros(Ny1,Nx1); % Old vy-velocity m/s
% P-nodes
RHO=zeros(Ny1,Nx1); % Density, kg/m^3
RHOCP=zeros(Ny1,Nx1); % Volumetric heat capacity, J/m^3/K
ALPHA=zeros(Ny1,Nx1); % Thermal expansion, 1/K
BETTA=zeros(Ny1,Nx1); % Compressibility, 1/Pa
HR=zeros(Ny1,Nx1); % Radioactive heating, W/m^3
HA=zeros(Ny1,Nx1); % Adiabatic heating, W/m^3
HS=zeros(Ny1,Nx1); % Shear heating, W/m^3
ETAP=zeros(Ny1,Nx1); % Viscosity, Pa*s
GGGP=zeros(Ny1,Nx1); % Shear modulus, Pa
EXX=zeros(Ny,Nx); % EPSILONxx, 1/s
SXX=zeros(Ny,Nx); % SIGMA'xx, 1/s
SXX0=zeros(Ny,Nx); % SIGMA0'xx, 1/s
tk1=zeros(Ny1,Nx1); % Old temperature, K
tk2=zeros(Ny1,Nx1); % New temperature, K
vxp=zeros(Ny1,Nx1); % Vx in pressure nodes, m/s
vyp=zeros(Ny1,Nx1); % Vy in pressure nodes, m/s
pr=zeros(Ny1,Nx1); % Pressure, Pa
pr0=zeros(Ny1,Nx1); % Old pressure, Pa

% Define markers
Nxm=(Nx-1)*4; % Marker grid resolution in horizontal direction
Nym=(Ny-1)*4; % Marker grid resolution in vertical direction
dxm=xsize/Nxm; % Marker grid step in horizontal direction,m
dym=ysize/Nym; % Marker grid step in vertical direction,m
marknum=Nxm*Nym; % Number of markers
xm=zeros(1,marknum); % Horizontal coordinates, m
ym=zeros(1,marknum); % Vertical coordinates, m
tm=zeros(1,marknum); % Material type
tkm=zeros(1,marknum); % Marker temperature, K
sxxm=zeros(1,marknum); % SIGMA'xx, Pa
sxym=zeros(1,marknum); % SIGMAxy, Pa
etavpm=zeros(1,marknum); % Visco-plastic viscosity, Pa
pr0m=zeros(1,marknum); % Marker pressure, Pa
vx0m=zeros(1,marknum); % Marker vx, m/s
vy0m=zeros(1,marknum); % Marker vy, m/s

% Define properties of materials: 
%        Block   Fault   Buffers
rhom    = [3000   3000   3000  ]; % Density, kg/m^3
etam    = [1e+23  1e+23  1e+23 ]; % Viscosity, Pa s
rhocpm  = [3.0e+6 3.0e+6 3.0e+6]; % Volumetric heat capacity, kg/m^3
alpham  = [3e-5   3e-5   3e-5  ]; % Thermal expansion, 1/K
bettam  = [1e-11  1e-11  1e-11 ]; % Compressibility, 1/Pa
km      = [3      3      3     ]; % Thermal conductivity, W/m/K
hrm     = [3e-8   3e-8   3e-8  ]; % Radiogenic heat production, W/m^3
gggm    = [1e+11  1e+11  1e+11 ]; % Shear Modulus, Pa
cohesm  = [3e+8   3e+7   3e+7  ]; % Cohesion, Pa
frictsm = [0.3    0.3    0.1   ]; % Static friction coefficient
frictdm = [0.3    0.1    0.3   ]; % Dynamic friction coefficient
Vc=1e-8; % Characteristic slip velocity, m/s

% Define marker coordinates, temperature and material type
m=1; % Marker counter
for jm=1:1:Nxm
    for im=1:1:Nym
        % Define marker coordinates
        xm(m)=dxm/2+(jm-1)*dxm+(rand-0.5)*dxm;
        ym(m)=dym/2+(im-1)*dym+(rand-0.5)*dym;
        % Marker properties
        % Block
        tm(m)=1; % Material type
        % Fault
        if(ym(m)>ysize*0.49 && ym(m)<ysize*0.51)
            tm(m)=2; % Material type$
            % Buffer zones 
            if(xm(m)<xsize*0.2 || xm(m)>xsize*0.8)
                tm(m)=3; % Material type
            end
        end
        % Initial temperature
        tkm(m)=273;
        % Initial visco-plastic viscosity
        etavpm(m)=etam(tm(m));
        % Update marker counter
        m=m+1;
    end
end
            
            
% Define global matrixes 
% Mechanical solution: L(), R()
N=Nx1*Ny1*3; % Global number of unknowns
L=sparse(N,N); % Matrix of coefficients (left part)
R=zeros(N,1); % Vector of right parts
% Thermal solution: LT(), RT()
N=Nx1*Ny1; % Global number of unknowns
LT=sparse(N,N); % Matrix of coefficients (left part)
RT=zeros(N,1); % Vector of right parts

% Mechanical boundary conditions: free slip=-1; No Slip=1
bcleft=1;
bcright=1;
bctop=1;
bcbottom=1;
% Top/Bottom velocities
strainrate=1e-13; % Shortening strain rate
vxtop=strainrate*xsize/2;
vxbottom=-strainrate*xsize/2;
% Pressure BC
PCONF=1e+8;

% Thermal boundary conditions: insulation at all boundaries

% Timestepping
dtelastic=1e+7; % Maximal computational timestep, s
dtmin=1e-6; % Minimal computational timestep, s
dt=dtelastic; % Current computational timestep, s
dtkoefv=1.1; % Koefficient to decrese dt for P,T,vx,vy,SIGMA limits 
dtkoef=2; % Koefficient to decrese dt in case of no convergence
dtkoefup=1.1; % Koefficient to increase dt
dtstep=200; % Number of iterations before changing dt
dxymax=0.0001; % Max marker movement per time step, grid steps
vpratio=1/3; % Weight of averaged velocity for moving markers
DTmax=20; % Max temperature change per time step, K
DSmax=1e+9; % Max stress change per time step, Pa
DPmax=1e+9; % Max pressure change per time step, Pa
dsubgridt=0; % Subgrid temperature diffusion parameter
dsubgrids=0; % Subgrid stress diffusion parameter
dsubgridv=0; % Subgrid velocity diffusion parameter
timesum=0; % Time sum, s
etamin=1e-3; % Lower viscosity cut-off, Pa s
etamax=1e+23; % Upper viscosity cut-off, Pa s
nplast=100000; % Number of plastic iterations
visstep=10; % Periodicity of visualization
yerrmax=3e+3; % Tolerance level for yielding error, Pa
terrmax=1e-3; % Tolerance level for temperature error, K
perrmax=3e+4; % Tolerance level for pressure error, Pa
YERRNOD=zeros(1,nplast); % Yielding error of nodes
TERRNOD=zeros(1,nplast); % Temperature error of nodes
PERRNOD=zeros(1,nplast); % Pressure error of nodes
etawt=0; % Weight for old viscosity
nsteps=5000; % number of timesteps
timestep=1;
for timestep=timestep:1:nsteps

% Save old stresses
sxxm00=sxxm; 
sxym00=sxym;    
    
% Interpolate properties from markers to nodes
% Basic nodes
ETA0SUM=zeros(Ny,Nx);
ETASUM=zeros(Ny,Nx);
GGGSUM=zeros(Ny,Nx);
SXYSUM=zeros(Ny,Nx);
COHSUM=zeros(Ny,Nx);
FRISSUM=zeros(Ny,Nx);
FRIDSUM=zeros(Ny,Nx);
WTSUM=zeros(Ny,Nx);
% Vx-nodes
RHOXSUM=zeros(Ny1,Nx1);
KXSUM=zeros(Ny1,Nx1);
VXSUM=zeros(Ny1,Nx1);
WTXSUM=zeros(Ny1,Nx1);
% Vy-nodes
RHOYSUM=zeros(Ny1,Nx1);
KYSUM=zeros(Ny1,Nx1);
VYSUM=zeros(Ny1,Nx1);
WTYSUM=zeros(Ny1,Nx1);
% P-Nodes
GGGPSUM=zeros(Ny1,Nx1);
SXXSUM=zeros(Ny1,Nx1);
RHOSUM=zeros(Ny1,Nx1);
RHOCPSUM=zeros(Ny1,Nx1);
ALPHASUM=zeros(Ny1,Nx1);
BETTASUM=zeros(Ny1,Nx1);
HRSUM=zeros(Ny1,Nx1);
TKSUM=zeros(Ny1,Nx1);
PRSUM=zeros(Ny1,Nx1);
WTPSUM=zeros(Ny1,Nx1);

for m=1:1:marknum
    
    % Interpolation to basic nodes
    % Define i,j indexes for the upper left node
    j=fix((xm(m)-x(1))/dx)+1;
    i=fix((ym(m)-y(1))/dy)+1;
    if(j<1)
        j=1;
    elseif(j>Nx-1)
        j=Nx-1;
    end
    if(i<1)
        i=1;
    elseif(i>Ny-1)
        i=Ny-1;
    end
    % Compute distances
    dxmj=xm(m)-x(j);
    dymi=ym(m)-y(i);
    % Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy);
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy);
    wtmi1j1=(dxmj/dx)*(dymi/dy);
    % Update properties
    % i,j Node
    ETA0SUM(i,j)=ETA0SUM(i,j)+etam(tm(m))*wtmij;
    ETASUM(i,j)=ETASUM(i,j)+etavpm(m)*wtmij;
    GGGSUM(i,j)=GGGSUM(i,j)+1/gggm(tm(m))*wtmij;
    SXYSUM(i,j)=SXYSUM(i,j)+sxym(m)*wtmij;
    COHSUM(i,j)=COHSUM(i,j)+cohesm(tm(m))*wtmij;
    FRISSUM(i,j)=FRISSUM(i,j)+frictsm(tm(m))*wtmij;
    FRIDSUM(i,j)=FRIDSUM(i,j)+frictdm(tm(m))*wtmij;
    WTSUM(i,j)=WTSUM(i,j)+wtmij;
    % i+1,j Node
    ETA0SUM(i+1,j)=ETA0SUM(i+1,j)+etam(tm(m))*wtmi1j;
    ETASUM(i+1,j)=ETASUM(i+1,j)+etavpm(m)*wtmi1j;
    GGGSUM(i+1,j)=GGGSUM(i+1,j)+1/gggm(tm(m))*wtmi1j;
    SXYSUM(i+1,j)=SXYSUM(i+1,j)+sxym(m)*wtmi1j;
    COHSUM(i+1,j)=COHSUM(i+1,j)+cohesm(tm(m))*wtmi1j;
    FRISSUM(i+1,j)=FRISSUM(i+1,j)+frictsm(tm(m))*wtmi1j;
    FRIDSUM(i+1,j)=FRIDSUM(i+1,j)+frictdm(tm(m))*wtmi1j;
    WTSUM(i+1,j)=WTSUM(i+1,j)+wtmi1j;
    % i,j+1 Node
    ETA0SUM(i,j+1)=ETA0SUM(i,j+1)+etam(tm(m))*wtmij1;
    ETASUM(i,j+1)=ETASUM(i,j+1)+etavpm(m)*wtmij1;
    GGGSUM(i,j+1)=GGGSUM(i,j+1)+1/gggm(tm(m))*wtmij1;
    SXYSUM(i,j+1)=SXYSUM(i,j+1)+sxym(m)*wtmij1;
    COHSUM(i,j+1)=COHSUM(i,j+1)+cohesm(tm(m))*wtmij1;
    FRISSUM(i,j+1)=FRISSUM(i,j+1)+frictsm(tm(m))*wtmij1;
    FRIDSUM(i,j+1)=FRIDSUM(i,j+1)+frictdm(tm(m))*wtmij1;
    WTSUM(i,j+1)=WTSUM(i,j+1)+wtmij1;
    % i+1,j+1 Node
    ETA0SUM(i+1,j+1)=ETA0SUM(i+1,j+1)+etam(tm(m))*wtmi1j1;
    ETASUM(i+1,j+1)=ETASUM(i+1,j+1)+etavpm(m)*wtmi1j1;
    GGGSUM(i+1,j+1)=GGGSUM(i+1,j+1)+1/gggm(tm(m))*wtmi1j1;
    SXYSUM(i+1,j+1)=SXYSUM(i+1,j+1)+sxym(m)*wtmi1j1;
    COHSUM(i+1,j+1)=COHSUM(i+1,j+1)+cohesm(tm(m))*wtmi1j1;
    FRISSUM(i+1,j+1)=FRISSUM(i+1,j+1)+frictsm(tm(m))*wtmi1j1;
    FRIDSUM(i+1,j+1)=FRIDSUM(i+1,j+1)+frictdm(tm(m))*wtmi1j1;
    WTSUM(i+1,j+1)=WTSUM(i+1,j+1)+wtmi1j1;    
    
    
    % Interpolation to vx-nodes
    % Define i,j indexes for the upper left node
    j=fix((xm(m)-xvx(1))/dx)+1;
    i=fix((ym(m)-yvx(1))/dy)+1;
    if(j<1)
        j=1;
    elseif(j>Nx-1)
        j=Nx-1;
    end
    if(i<1)
        i=1;
    elseif(i>Ny)
        i=Ny;
    end
    % Compute distances
    dxmj=xm(m)-xvx(j);
    dymi=ym(m)-yvx(i);
    % Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy);
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy);
    wtmi1j1=(dxmj/dx)*(dymi/dy);
    % Update properties
    % i,j Node
    RHOXSUM(i,j)=RHOXSUM(i,j)+rhom(tm(m))*wtmij;
    KXSUM(i,j)=KXSUM(i,j)+km(tm(m))*wtmij;
    VXSUM(i,j)=VXSUM(i,j)+vx0m(m)*rhom(tm(m))*wtmij;
    WTXSUM(i,j)=WTXSUM(i,j)+wtmij;
    % i+1,j Node
    RHOXSUM(i+1,j)=RHOXSUM(i+1,j)+rhom(tm(m))*wtmi1j;
    KXSUM(i+1,j)=KXSUM(i+1,j)+km(tm(m))*wtmi1j;
    VXSUM(i+1,j)=VXSUM(i+1,j)+vx0m(m)*rhom(tm(m))*wtmi1j;
    WTXSUM(i+1,j)=WTXSUM(i+1,j)+wtmi1j;
    % i,j+1 Node
    RHOXSUM(i,j+1)=RHOXSUM(i,j+1)+rhom(tm(m))*wtmij1;
    KXSUM(i,j+1)=KXSUM(i,j+1)+km(tm(m))*wtmij1;
    VXSUM(i,j+1)=VXSUM(i,j+1)+vx0m(m)*rhom(tm(m))*wtmij1;
    WTXSUM(i,j+1)=WTXSUM(i,j+1)+wtmij1;
    % i+1,j+1 Node
    RHOXSUM(i+1,j+1)=RHOXSUM(i+1,j+1)+rhom(tm(m))*wtmi1j1;
    KXSUM(i+1,j+1)=KXSUM(i+1,j+1)+km(tm(m))*wtmi1j1;
    VXSUM(i+1,j+1)=VXSUM(i+1,j+1)+vx0m(m)*rhom(tm(m))*wtmi1j1;
    WTXSUM(i+1,j+1)=WTXSUM(i+1,j+1)+wtmi1j1;

    % Interpolation to vy-nodes
    % Define i,j indexes for the upper left node
    j=fix((xm(m)-xvy(1))/dx)+1;
    i=fix((ym(m)-yvy(1))/dy)+1;
    if(j<1)
        j=1;
    elseif(j>Nx)
        j=Nx;
    end
    if(i<1)
        i=1;
    elseif(i>Ny-1)
        i=Ny-1;
    end
    % Compute distances
    dxmj=xm(m)-xvy(j);
    dymi=ym(m)-yvy(i);
    % Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy);
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy);
    wtmi1j1=(dxmj/dx)*(dymi/dy);
    % Update properties
    % i,j Node
    RHOYSUM(i,j)=RHOYSUM(i,j)+rhom(tm(m))*wtmij;
    KYSUM(i,j)=KYSUM(i,j)+km(tm(m))*wtmij;
    VYSUM(i,j)=VYSUM(i,j)+vy0m(m)*rhom(tm(m))*wtmij;
    WTYSUM(i,j)=WTYSUM(i,j)+wtmij;
    % i+1,j Node
    RHOYSUM(i+1,j)=RHOYSUM(i+1,j)+rhom(tm(m))*wtmi1j;
    KYSUM(i+1,j)=KYSUM(i+1,j)+km(tm(m))*wtmi1j;
    VYSUM(i+1,j)=VYSUM(i+1,j)+vy0m(m)*rhom(tm(m))*wtmi1j;
    WTYSUM(i+1,j)=WTYSUM(i+1,j)+wtmi1j;
    % i,j+1 Node
    RHOYSUM(i,j+1)=RHOYSUM(i,j+1)+rhom(tm(m))*wtmij1;
    KYSUM(i,j+1)=KYSUM(i,j+1)+km(tm(m))*wtmij1;
    VYSUM(i,j+1)=VYSUM(i,j+1)+vy0m(m)*rhom(tm(m))*wtmij1;
    WTYSUM(i,j+1)=WTYSUM(i,j+1)+wtmij1;
    % i+1,j+1 Node
    RHOYSUM(i+1,j+1)=RHOYSUM(i+1,j+1)+rhom(tm(m))*wtmi1j1;
    KYSUM(i+1,j+1)=KYSUM(i+1,j+1)+km(tm(m))*wtmi1j1;
    VYSUM(i+1,j+1)=VYSUM(i+1,j+1)+vy0m(m)*rhom(tm(m))*wtmi1j1;
    WTYSUM(i+1,j+1)=WTYSUM(i+1,j+1)+wtmi1j1;

    % Interpolation to P-nodes
    % Define i,j indexes for the upper left node
    j=fix((xm(m)-xp(1))/dx)+1;
    i=fix((ym(m)-yp(1))/dy)+1;
    if(j<1)
        j=1;
    elseif(j>Nx)
        j=Nx;
    end
    if(i<1)
        i=1;
    elseif(i>Ny)
        i=Ny;
    end
    % Compute distances
    dxmj=xm(m)-xp(j);
    dymi=ym(m)-yp(i);
    % Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy);
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy);
    wtmi1j1=(dxmj/dx)*(dymi/dy);
    % Update properties
    % i,j Node
    GGGPSUM(i,j)=GGGPSUM(i,j)+1/gggm(tm(m))*wtmij;
    SXXSUM(i,j)=SXXSUM(i,j)+sxxm(m)*wtmij;
    RHOSUM(i,j)=RHOSUM(i,j)+rhom(tm(m))*wtmij;
    RHOCPSUM(i,j)=RHOCPSUM(i,j)+rhocpm(tm(m))*wtmij;
    ALPHASUM(i,j)=ALPHASUM(i,j)+alpham(tm(m))*wtmij;
    BETTASUM(i,j)=BETTASUM(i,j)+bettam(tm(m))*wtmij;
    HRSUM(i,j)=HRSUM(i,j)+hrm(tm(m))*wtmij;
    TKSUM(i,j)=TKSUM(i,j)+tkm(m)*rhocpm(tm(m))*wtmij;
    PRSUM(i,j)=PRSUM(i,j)+pr0m(m)*wtmij;
    WTPSUM(i,j)=WTPSUM(i,j)+wtmij;
    % i+1,j Node
    GGGPSUM(i+1,j)=GGGPSUM(i+1,j)+1/gggm(tm(m))*wtmi1j;
    SXXSUM(i+1,j)=SXXSUM(i+1,j)+sxxm(m)*wtmi1j;
    RHOSUM(i+1,j)=RHOSUM(i+1,j)+rhom(tm(m))*wtmi1j;
    RHOCPSUM(i+1,j)=RHOCPSUM(i+1,j)+rhocpm(tm(m))*wtmi1j;
    ALPHASUM(i+1,j)=ALPHASUM(i+1,j)+alpham(tm(m))*wtmi1j;
    BETTASUM(i+1,j)=BETTASUM(i+1,j)+bettam(tm(m))*wtmi1j;
    HRSUM(i+1,j)=HRSUM(i+1,j)+hrm(tm(m))*wtmi1j;
    TKSUM(i+1,j)=TKSUM(i+1,j)+tkm(m)*rhocpm(tm(m))*wtmi1j;
    PRSUM(i+1,j)=PRSUM(i+1,j)+pr0m(m)*wtmi1j;
    WTPSUM(i+1,j)=WTPSUM(i+1,j)+wtmi1j;
    % i,j+1 Node
    GGGPSUM(i,j+1)=GGGPSUM(i,j+1)+1/gggm(tm(m))*wtmij1;
    SXXSUM(i,j+1)=SXXSUM(i,j+1)+sxxm(m)*wtmij1;
    RHOSUM(i,j+1)=RHOSUM(i,j+1)+rhom(tm(m))*wtmij1;
    RHOCPSUM(i,j+1)=RHOCPSUM(i,j+1)+rhocpm(tm(m))*wtmij1;
    ALPHASUM(i,j+1)=ALPHASUM(i,j+1)+alpham(tm(m))*wtmij1;
    BETTASUM(i,j+1)=BETTASUM(i,j+1)+bettam(tm(m))*wtmij1;
    HRSUM(i,j+1)=HRSUM(i,j+1)+hrm(tm(m))*wtmij1;
    TKSUM(i,j+1)=TKSUM(i,j+1)+tkm(m)*rhocpm(tm(m))*wtmij1;
    PRSUM(i,j+1)=PRSUM(i,j+1)+pr0m(m)*wtmij1;
    WTPSUM(i,j+1)=WTPSUM(i,j+1)+wtmij1;
    % i+1,j+1 Node
    GGGPSUM(i+1,j+1)=GGGPSUM(i+1,j+1)+1/gggm(tm(m))*wtmi1j1;
    SXXSUM(i+1,j+1)=SXXSUM(i+1,j+1)+sxxm(m)*wtmi1j1;
    RHOSUM(i+1,j+1)=RHOSUM(i+1,j+1)+rhom(tm(m))*wtmi1j1;
    RHOCPSUM(i+1,j+1)=RHOCPSUM(i+1,j+1)+rhocpm(tm(m))*wtmi1j1;
    ALPHASUM(i+1,j+1)=ALPHASUM(i+1,j+1)+alpham(tm(m))*wtmi1j1;
    BETTASUM(i+1,j+1)=BETTASUM(i+1,j+1)+bettam(tm(m))*wtmi1j1;
    HRSUM(i+1,j+1)=HRSUM(i+1,j+1)+hrm(tm(m))*wtmi1j1;
    TKSUM(i+1,j+1)=TKSUM(i+1,j+1)+tkm(m)*rhocpm(tm(m))*wtmi1j1;
    PRSUM(i+1,j+1)=PRSUM(i+1,j+1)+pr0m(m)*wtmi1j1;
    WTPSUM(i+1,j+1)=WTPSUM(i+1,j+1)+wtmi1j1;
end
% Compute physical properties
% Basic nodes
% YNY=zeros(Ny,Nx);
for j=1:1:Nx
    for i=1:1:Ny
        if(WTSUM(i,j)>0)
            ETA0(i,j)=ETA0SUM(i,j)/WTSUM(i,j);
%             ETA(i,j)=ETASUM(i,j)/WTSUM(i,j);
%             if(ETA(i,j)<ETA0(i,j))
%                 YNY(i,j)=1;
%             end
            GGG(i,j)=1/(GGGSUM(i,j)/WTSUM(i,j));
%             SXY0(i,j)=SXYSUM(i,j)/WTSUM(i,j);
            COH(i,j)=COHSUM(i,j)/WTSUM(i,j);
            FRIS(i,j)=FRISSUM(i,j)/WTSUM(i,j);
            FRID(i,j)=FRIDSUM(i,j)/WTSUM(i,j);
        end
    end
end
% Vx-nodes
for j=1:1:Nx1
    for i=1:1:Ny1
        if(WTXSUM(i,j)>0)
            RHOX(i,j)=RHOXSUM(i,j)/WTXSUM(i,j);
            KX(i,j)=KXSUM(i,j)/WTXSUM(i,j);
%             vx0(i,j)=VXSUM(i,j)/RHOXSUM(i,j);
        end
    end
end
% Vy-nodes
for j=1:1:Nx1
    for i=1:1:Ny1
        if(WTYSUM(i,j)>0)
            RHOY(i,j)=RHOYSUM(i,j)/WTYSUM(i,j);
            KY(i,j)=KYSUM(i,j)/WTYSUM(i,j);
%             vy0(i,j)=VYSUM(i,j)/RHOYSUM(i,j);
        end
    end
end
% P-nodes
for j=1:1:Nx1
    for i=1:1:Ny1
        if(WTPSUM(i,j)>0)
            GGGP(i,j)=1/(GGGPSUM(i,j)/WTPSUM(i,j));
%             SXX0(i,j)=SXXSUM(i,j)/WTPSUM(i,j);
            RHO(i,j)=RHOSUM(i,j)/WTPSUM(i,j);
            RHOCP(i,j)=RHOCPSUM(i,j)/WTPSUM(i,j);
            ALPHA(i,j)=ALPHASUM(i,j)/WTPSUM(i,j);
            BETTA(i,j)=BETTASUM(i,j)/WTPSUM(i,j);
            HR(i,j)=HRSUM(i,j)/WTPSUM(i,j);
            tk1(i,j)=TKSUM(i,j)/RHOCPSUM(i,j);
%             pr0(i,j)=PRSUM(i,j)/WTPSUM(i,j);
        end
    end
end
% Applying thermal boundary conditions for interpolated temperature
% Upper boundary 
tk1(1,2:Nx)=tk1(2,2:Nx); % Insulating boundary
% Lower boundary 
tk1(Ny1,2:Nx)=tk1(Ny,2:Nx); % Insulating boundary
% Left boundary
tk1(:,1)=tk1(:,2); % Insulating boundary
% Right boundary
tk1(:,Nx1)=tk1(:,Nx); % Insulating boundary
% Applying vx-velocity boundary conditions for interpolated velocity
% Upper boundary 
vx0(1,2:Nx-1)=2*vxtop-bctop*vx0(2,2:Nx-1); % prescribed velocity
% Lower boundary 
vx0(Ny1,2:Nx-1)=2*vxbottom-bcbottom*vx0(Ny,2:Nx-1); % prescribed velocity
% Left boundary
vx0(:,1)=vx0(:,2); % dVx/dx=0
% Right boundary
vx0(:,Nx)=vx0(:,Nx-1); % dVx/dx=0
% Applying vy-velocity boundary conditions for interpolated velocity
% Left boundary
vy0(2:Ny-1,1)=-bcleft*vy0(2:Ny-1,2); % no slip
% Right boundary
vy0(2:Ny-1,Nx1)=-bcright*vy0(2:Ny-1,Nx); % no slip
% Upper boundary 
vy0(1,:)=0; % no slip
% Lower boundary 
vy0(Ny,:)=0; % no slip
% Applying Pressure boundary conditions for interpolated pressure
dpr0=pr0(2,2)-PCONF;
pr0=pr0-dpr0;
pr0m=pr0m-dpr0;

% Try to increase computational Timestep
dt=min(dt*dtkoefup,dtelastic)

% Set initial viscoplastic viscosity
if(timestep==1)
    ETA=ETA0;
end


% Save initial viscoplastic viscosity
ETA00=ETA;
% Save initial yielding nodes
YNY00=YNY;
% Set initial new temperature
tk2=tk1;

% Start thermomechanical iterations on Nodes
icount=0; % number of iterations without changing timestep
for iplast=1:1:nplast
% Set lower limit for viscosity in Navier-Stokes equations
etaminNS=max(etamin, dt*min(min(GGG))*1e-9);
% Limit viscosity in the basic nodes
for i=1:1:Ny
    for j=1:1:Nx
        if(ETA(i,j)>ETA0(i,j))
            ETA(i,j)=ETA0(i,j);
        end
        if(ETA(i,j)<etaminNS)
            ETA(i,j)=etaminNS;
        end
    end
end
% Introducing scaled pressure
pscale=etaminNS/dx;
% Recompute viscosity at pressure nodes
for i=2:1:Ny
    for j=2:1:Nx
        ETAP(i,j)=1/((1/ETA(i-1,j-1)+1/ETA(i,j-1)+1/ETA(i-1,j)+1/ETA(i,j))/4);
    end
end
% Composing global matrixes L(), R() for Navier-Stokes and continuity equations
for j=1:1:Nx1
    for i=1:1:Ny1
        % Define global indexes in algebraic space
        kvx=((j-1)*Ny1+i-1)*3+1; % Vx
        kvy=kvx+1; % Vy
        kpm=kvx+2; % P
        
        % Vx equation External points
        if(i==1 || i==Ny1 || j==1 || j==Nx || j==Nx1)
            % Boundary Condition 
            % Ghost unknowns 1*Vx=0
            if(j==Nx1)
                L(kvx,kvx)=1; % Left part
                R(kvx)=0; % Right part
            end
            % Left Boundary
            if(j==1)
                L(kvx,kvx)=1; % Left part
                L(kvx,kvx+Ny1*3)=-1; % Left part
                R(kvx)=0; % Right part
            end
            % Right Boundary
            if(j==Nx)
                L(kvx,kvx)=1; % Left part
                L(kvx,kvx-Ny1*3)=-1; % Left part
                R(kvx)=0; % Right part
            end
            % Top boundary
            if(i==1 && j>1 && j<Nx)
                L(kvx,kvx)=1; % Left part
                L(kvx,kvx+3)=bctop; % Left part
                R(kvx)=2*vxtop; % Right part
            end
            % Top boundary
            if(i==Ny1 && j>1 && j<Nx)
                L(kvx,kvx)=1; % Left part
                L(kvx,kvx-3)=bcbottom; % Left part
                R(kvx)=2*vxbottom; % Right part
            end
        else
        % Internal points: x-Navier-Stokes eq.
        %            Vx2
        %             |
        %        Vy1  |  Vy3
        %             |
        %     Vx1-P1-Vx3-P2-Vx5
        %             |
        %        Vy2  |  Vy4
        %             |
        %            Vx4
        %
        % Computational viscosity
        ETA1=ETA(i-1,j)*GGG(i-1,j)*dt/(GGG(i-1,j)*dt+ETA(i-1,j));
        ETA2=ETA(i,j)*GGG(i,j)*dt/(GGG(i,j)*dt+ETA(i,j));
        ETAP1=ETAP(i,j)*GGGP(i,j)*dt/(GGGP(i,j)*dt+ETAP(i,j));
        ETAP2=ETAP(i,j+1)*GGGP(i,j+1)*dt/(GGGP(i,j+1)*dt+ETAP(i,j+1));
        % Old stresses
        SXY1=SXY0(i-1,j)*ETA(i-1,j)/(GGG(i-1,j)*dt+ETA(i-1,j));
        SXY2=SXY0(i,j)*ETA(i,j)/(GGG(i,j)*dt+ETA(i,j));
        SXX1=SXX0(i,j)*ETAP(i,j)/(GGGP(i,j)*dt+ETAP(i,j));
        SXX2=SXX0(i,j+1)*ETAP(i,j+1)/(GGGP(i,j+1)*dt+ETAP(i,j+1));
        % Density gradients
        dRHOdx=(RHOX(i,j+1)-RHOX(i,j-1))/2/dx;
        dRHOdy=(RHOX(i+1,j)-RHOX(i-1,j))/2/dy;
        % Left part
        L(kvx,kvx-Ny1*3)=ETAP1/dx^2; % Vx1
        L(kvx,kvx-3)=ETA1/dy^2; % Vx2
        L(kvx,kvx)=-(ETAP1+ETAP2)/dx^2-...
                      (ETA1+ETA2)/dy^2-...
                      dRHOdx*gx*dt-RHOX(i,j)/dt; % Vx3
        L(kvx,kvx+3)=ETA2/dy^2; % Vx4
        L(kvx,kvx+Ny1*3)=ETAP2/dx^2; % Vx5
        L(kvx,kvy)=ETAP1/dx/dy-ETA2/dx/dy-dRHOdy*gx*dt/4;  % Vy2
        L(kvx,kvy+Ny1*3)=-ETAP2/dx/dy+ETA2/dx/dy-dRHOdy*gx*dt/4;  % Vy4
        L(kvx,kvy-3)=-ETAP1/dx/dy+ETA1/dx/dy-dRHOdy*gx*dt/4;  % Vy1
        L(kvx,kvy+Ny1*3-3)=ETAP2/dx/dy-ETA1/dx/dy-dRHOdy*gx*dt/4;  % Vy3
        L(kvx,kpm)=pscale/dx; % P1
        L(kvx,kpm+Ny1*3)=-pscale/dx; % P2
        % Right part
        R(kvx)=-RHOX(i,j)*(gx+vx0(i,j)/dt)-(SXY2-SXY1)/dy-(SXX2-SXX1)/dx;
        end
        
        % Vy equation External points
        if(j==1 || j==Nx1 || i==1 || i==Ny || i==Ny1)
            % Boundary Condition
            % Ghost unknowns 1*Vx=0
            if(i==Ny1)
                L(kvy,kvy)=1; % Left part
                R(kvy)=0; % Right part
            end
            % Top boundary
            if(i==1)
                L(kvy,kvy)=1; % Left part
                R(kvy)=0; % Right part
            end
            % Bottom boundary
            if(i==Ny)
                L(kvy,kvy)=1; % Left part
                R(kvy)=0; % Right part
            end
            % Left boundary
            if(j==1 && i>1 && i<Ny)
                L(kvy,kvy)=1; % Left part
                L(kvy,kvy+3*Ny1)=bcleft; % Left part
                R(kvy)=0; % Right part
            end
            % Right boundary
            if(j==Nx1 && i>1 && i<Ny)
                L(kvy,kvy)=1; % Left part
                L(kvy,kvy-3*Ny1)=bcright; % Left part
                R(kvy)=0; % Right part
            end
        else
        % Internal points: y-Navier-Stokes eq.
        %            Vy2
        %             |
        %         Vx1 P1 Vx3
        %             |
        %     Vy1----Vy3----Vy5
        %             |
        %         Vx2 P2 Vx4
        %             |
        %            Vy4
        %
        % Computational viscosity
        ETA1=ETA(i,j-1)*GGG(i,j-1)*dt/(GGG(i,j-1)*dt+ETA(i,j-1));
        ETA2=ETA(i,j)*GGG(i,j)*dt/(GGG(i,j)*dt+ETA(i,j));
        ETAP1=ETAP(i,j)*GGGP(i,j)*dt/(GGGP(i,j)*dt+ETAP(i,j));
        ETAP2=ETAP(i+1,j)*GGGP(i+1,j)*dt/(GGGP(i+1,j)*dt+ETAP(i+1,j));
        % Old stresses
        SXY1=SXY0(i,j-1)*ETA(i,j-1)/(GGG(i,j-1)*dt+ETA(i,j-1));
        SXY2=SXY0(i,j)*ETA(i,j)/(GGG(i,j)*dt+ETA(i,j));
        SYY1=-SXX0(i,j)*ETAP(i,j)/(GGGP(i,j)*dt+ETAP(i,j));
        SYY2=-SXX0(i+1,j)*ETAP(i+1,j)/(GGGP(i+1,j)*dt+ETAP(i+1,j));
        % Density gradients
        dRHOdx=(RHOY(i,j+1)-RHOY(i,j-1))/2/dx;
        dRHOdy=(RHOY(i+1,j)-RHOY(i-1,j))/2/dy;
        % Left part
        L(kvy,kvy-Ny1*3)=ETA1/dx^2; % Vy1
        L(kvy,kvy-3)=ETAP1/dy^2; % Vy2
        L(kvy,kvy)=-(ETAP1+ETAP2)/dy^2-...
                      (ETA1+ETA2)/dx^2-...
                      dRHOdy*gy*dt-RHOY(i,j)/dt; % Vy3
        L(kvy,kvy+3)=ETAP2/dy^2; % Vy4
        L(kvy,kvy+Ny1*3)=ETA2/dx^2; % Vy5
        L(kvy,kvx)=ETAP1/dx/dy-ETA2/dx/dy-dRHOdx*gy*dt/4; %Vx3
        L(kvy,kvx+3)=-ETAP2/dx/dy+ETA2/dx/dy-dRHOdx*gy*dt/4; %Vx4
        L(kvy,kvx-Ny1*3)=-ETAP1/dx/dy+ETA1/dx/dy-dRHOdx*gy*dt/4; %Vx1
        L(kvy,kvx+3-Ny1*3)=ETAP2/dx/dy-ETA1/dx/dy-dRHOdx*gy*dt/4; %Vx2
        L(kvy,kpm)=pscale/dy; % P1
        L(kvy,kpm+3)=-pscale/dy; % P2
        
        % Right part
        R(kvy)=-RHOY(i,j)*(gy+vy0(i,j)/dt)-(SXY2-SXY1)/dx-(SYY2-SYY1)/dy;
        end
        
        % P equation External points
        if(i==1 || j==1 || i==Ny1 || j==Nx1)
            % Boundary Condition
            % 1*P=0
            L(kpm,kpm)=1; % Left part
            R(kpm)=0; % Right part
        else
        % Internal points: compressible continuity eq.
        % dVx/dx+dVy/dy+BETTA*DP/Dt=ALPHA*DT/Dt
        %            Vy1
        %             |
        %        Vx1--P--Vx2
        %             |
        %            Vy2
        %
        % Left part
        L(kpm,kvx-Ny1*3)=-1/dx; % Vx1
        L(kpm,kvx)=1/dx; % Vx2
        L(kpm,kvy-3)=-1/dy; % Vy1
        L(kpm,kvy)=1/dy; % Vy2
        % Incompressible pressure buildup at the first timestep
%         if(timestep>1)
            L(kpm,kpm)=pscale*BETTA(i,j)/dt; % P
%         end
        % Right part
        % Incompressible pressure buildup at the first timestep
%         if(timestep>1)
            R(kpm)=BETTA(i,j)*pr0(i,j)/dt+ALPHA(i,j)*(tk2(i,j)-tk1(i,j))/dt;
%         else
%             R(kpm)=0;
%         end
            
        end
        
    end
end

% 4) Solving matrixes, reloading solution
S=L\R; % Obtaining algebraic vector of solutions S()
% Reload solutions S() to vx(), vy(), p()
% Going through all grid points
for j=1:1:Nx1
    for i=1:1:Ny1
        % Define global indexes in algebraic space
        kvx=((j-1)*Ny1+i-1)*3+1; % Vx
        kvy=kvx+1; % Vy
        kpm=kvx+2; % P
        % Reload solution
        vx(i,j)=S(kvx);
        vy(i,j)=S(kvy);
        pr(i,j)=S(kpm)*pscale;
    end
end
% if(timestep==1)
%     pr0=pr;
% end

% Compute DVX,DVY
DVX=vx-vx0;
DVY=vy-vy0;


% Compute Stress, stress change and strain rate components
% Compute EPSILONxy, SIGMAxy in basic nodes
EXY=zeros(Ny,Nx); % Strain rate EPSILONxy, 1/s
SXY=zeros(Ny,Nx); % Stress SIGMAxy, Pa
DSXY=zeros(Ny,Nx); % Stress change SIGMAxy, Pa
for j=1:1:Nx
    for i=1:1:Ny
        % EXY,SXY, DSXY
        EXY(i,j)=0.5*((vx(i+1,j)-vx(i,j))/dy+...
            (vy(i,j+1)-vy(i,j))/dx);
        SXY(i,j)=2*ETA(i,j)*EXY(i,j)*GGG(i,j)*dt/(GGG(i,j)*dt+ETA(i,j))+...
            SXY0(i,j)*ETA(i,j)/(GGG(i,j)*dt+ETA(i,j));
        DSXY(i,j)=SXY(i,j)-SXY0(i,j);
    end
end
% Compute EPSILONxx, SIGMA'xx in pressure nodes
EXX=zeros(Ny1,Nx1); % Strain rate EPSILONxx, 1/s
EII=zeros(Ny1,Nx1); % Second strain rate invariant, 1/s
SXX=zeros(Ny1,Nx1); % Stress SIGMA'xx, Pa
SII=zeros(Ny1,Nx1); % Second stress invariant, Pa
DSXX=zeros(Ny1,Nx1); % Stress change SIGMA'xx, Pa
DP=zeros(Ny1,Nx1); % Pressure change, Pa
for j=2:1:Nx
    for i=2:1:Ny
        % EXX
        EXX(i,j)=0.5*((vx(i,j)-vx(i,j-1))/dx-(vy(i,j)-vy(i-1,j))/dy);
        % SXX
        SXX(i,j)=2*ETAP(i,j)*EXX(i,j)*GGGP(i,j)*dt/(GGGP(i,j)*dt+ETAP(i,j))+...
            SXX0(i,j)*ETAP(i,j)/(GGGP(i,j)*dt+ETAP(i,j));
        DSXX(i,j)=SXX(i,j)-SXX0(i,j);
        % EII
        EII(i,j)=(EXX(i,j)^2+((EXY(i,j)+EXY(i-1,j)+...
            EXY(i,j-1)+EXY(i-1,j-1))/4)^2)^0.5;
        % SII
        SII(i,j)=(SXX(i,j)^2+((SXY(i,j)+SXY(i-1,j)+...
            SXY(i,j-1)+SXY(i-1,j-1))/4)^2)^0.5;
        % DP
        DP(i,j)=pr(i,j)-pr0(i,j);
    end
end


% Compute shear heating HS in Temperature/Pressure nodes
HS=zeros(Ny1,Nx1); % Adiabatic heating, W/m^3
for j=2:1:Nx
    for i=2:1:Ny
        % Average SXY*EXY
        SXYEXY=(SXY(i,j)^2/ETA(i,j)+SXY(i-1,j)^2/ETA(i-1,j)+...
            SXY(i,j-1)^2/ETA(i,j-1)+SXY(i-1,j-1)^2/ETA(i-1,j-1))/4;
        % HS
        HS(i,j)=SXX(i,j)^2/ETAP(i,j)+SXYEXY;
    end
end

% Compute adiabatic heating HA in Temperature/Pressure nodes
HA=zeros(Ny1,Nx1); % Shear heating, W/m^3
for j=2:1:Nx
    for i=2:1:Ny
        % Average vy, vx
        VXP=(vx(i,j)+vx(i,j-1))/2;
        VYP=(vy(i,j)+vy(i-1,j))/2;
        % HA
        HA(i,j)=tk1(i,j)*ALPHA(i,j)*(pr(i,j)-pr0(i,j))/dt;
    end
end

% Composing global matrixes LT(), RT()
% Going through all points of the 2D grid and
% composing respective equations
for j=1:1:Nx1
    for i=1:1:Ny1
        % Define global index in algebraic space
        gk=(j-1)*Ny1+i;
        % External points
        if(i==1 || i==Ny1 || j==1 || j==Nx1)
            % Boundary Condition
            % Top BC: T=273
            if(i==1 && j>1 && j<Nx1)
                LT(gk,gk)=1; % Left part
                LT(gk,gk+1)=-1; % Left part
                RT(gk)=0; % Right part
            end
            % Bottom BC: T=1500
            if(i==Ny1 && j>1 && j<Nx1)
                LT(gk,gk)=1; % Left part
                LT(gk,gk-1)=-1; % Left part
                RT(gk)=0; % Right part
            end
            % Left BC: dT/dx=0
            if(j==1)
                LT(gk,gk)=1; % Left part
                LT(gk,gk+Ny1)=-1; % Left part
                RT(gk)=0; % Right part
            end
            % Right BC: dT/dx=0
            if(j==Nx1)
                LT(gk,gk)=1; % Left part
                LT(gk,gk-Ny1)=-1; % Left part
                RT(gk)=0; % Right part
            end
        else
        % Internal points: Temperature eq.
        % RHO*CP*dT/dt=-dqx/dx-dqy/dy+Hr+Hs+Ha
        %          Tdt2
        %           |
        %          Ky1
        %           |
        %Tdt1-Kx1-T03,Tdt3-Kx2-Tdt5
        %           |
        %          Ky2
        %           |
        %          Tdt4
        %
        % Left part
        Kx1=KX(i,j-1); 
        Kx2=KX(i,j); 
        Ky1=KY(i-1,j); 
        Ky2=KY(i,j); 
        LT(gk,gk-Ny1)=-Kx1/dx^2; % T1
        LT(gk,gk-1)=-Ky1/dy^2; % FI2
        LT(gk,gk)=RHOCP(i,j)/dt+(Kx1+Kx2)/dx^2+(Ky1+Ky2)/dy^2; % FI3
        LT(gk,gk+1)=-Ky2/dy^2; % FI4
        LT(gk,gk+Ny1)=-Kx2/dx^2; % FI5
        % Right part
        RT(gk)=RHOCP(i,j)/dt*tk1(i,j)+HR(i,j)+HA(i,j)+HS(i,j);
        end
    end
end

% Solving matrixes
ST=LT\RT; % Obtaining algebraic vector of solutions ST()

% Reload solutions ST() to geometrical array Tdt()
% Going through all grid points
for j=1:1:Nx1
    for i=1:1:Ny1
        % Compute global index
        gk=(j-1)*Ny1+i;
        % Reload solution
        tk2(i,j)=ST(gk);
    end
end
% Compute DT
DT=tk2-tk1;

% Nodal adjusment
% External P-nodes: symmetry
pr(:,1)=pr(:,2);
pr(:,Nx1)=pr(:,Nx);
pr(1,:)=pr(2,:);
pr(Ny1,:)=pr(Ny,:);
pr0(:,1)=pr0(:,2);
pr0(:,Nx1)=pr0(:,Nx);
pr0(1,:)=pr0(2,:);
pr0(Ny1,:)=pr0(Ny,:);
SXX(:,1)=SXX(:,2);
SXX(:,Nx1)=SXX(:,Nx);
SXX(1,:)=SXX(2,:);
SXX(Ny1,:)=SXX(Ny,:);
SXX0(:,1)=SXX0(:,2);
SXX0(:,Nx1)=SXX0(:,Nx);
SXX0(1,:)=SXX0(2,:);
SXX0(Ny1,:)=SXX0(Ny,:);
% Update viscosity for yielding
% Basic nodes
ETA5=ETA0;
YNY5=zeros(Ny,Nx);
DSY=zeros(Ny,Nx);
ynpl=0;
ddd=0;
for i=1:1:Ny
  for j=1:1:Nx
    % Compute second stress invariant
    % SXX, pt are averaged from four surrounding pressure nodes
    SIIB=(SXY(i,j)^2+((SXX(i,j)+SXX(i+1,j)+SXX(i,j+1)+SXX(i+1,j+1))/4)^2)^0.5;
    % Compute second invariant for a purely elastic stress build-up
    siiel=SIIB*(GGG(i,j)*dt+ETA(i,j))/ETA(i,j);
    % Compute pressure
    prB=(pr(i,j)+pr(i+1,j)+pr(i,j+1)+pr(i+1,j+1))/4;
    % Compute plastic slip rate
    V=YNY(i,j)*SIIB*(1/ETA(i,j)-1/ETA0(i,j))*dx;
    % Compute dynamic yielding stress
    syield=COH(i,j)+(FRID(i,j)+(FRIS(i,j)-FRID(i,j))/(1+V/Vc))*prB;
    syield=max(syield,0); % Non-negative strength requirement
    syieldit=0;
    etaplit=0;
    for ri=1:1:3
        etapl=ETA0(i,j);
        if(syield<siiel)
            etapl=max(etamin,dt*GGG(i,j)*syield/(siiel-syield));
            if(etapl<ETA0(i,j))
                V=syield*(1/etapl-1/ETA0(i,j))*dx;
%                 V=SIIB*(1/etapl-1/ETA0(i,j))*dx;
%                 V=(syield*SIIB)^0.5*(1/etapl-1/ETA0(i,j))*dx;
            else
                V=0;
                etapl=ETA0(i,j);
            end
        end
        syield=COH(i,j)+(FRID(i,j)+(FRIS(i,j)-FRID(i,j))/(1+V/Vc))*prB;
        syield=max(syield,0); % Non-negative strength requirement
        syieldit(ri)=syield;
        etaplit(ri)=etapl;
    end
%     if(YNY(i,j)>0)
%         figure(5);clf
%         subplot(2,1,1)
%         plot(log10(etaplit))
%         title(['j=',num2str(j),' i=',num2str(i)])
%         subplot(2,1,2)
%         plot(log10(syieldit))
%         pause(1)
%     end
    % Update error for old yielding nodes
    ynn=0;
    if(YNY(i,j)>0)
        ynn=1;
        DSY(i,j)=SIIB-syield;
        ddd=ddd+DSY(i,j)^2;
        ynpl=ynpl+1;
    end
    % Correcting viscosity for yielding
    if(syield<siiel)
        % New viscosity for the basic node
        etapl=dt*GGG(i,j)*syield/(siiel-syield);
        if(etapl<ETA0(i,j))
            % Recompute nodal visocity
            ETA5(i,j)=etapl^(1-etawt)*ETA(i,j)^etawt;
            % Mark yielding nodes
            YNY5(i,j)=1;
            % Apply viscosity cutoff values
            if(ETA5(i,j)<etamin)
                ETA5(i,j)=etamin;
            elseif(ETA5(i,j)>etamax)
                ETA5(i,j)=etamax;
            end
            % Update Error for new yielding nodes
            if(ynn==0)
                DSY(i,j)=SIIB-syield;
                ddd=ddd+DSY(i,j)^2;
                ynpl=ynpl+1;
            end
        else
            ETA5(i,j)=ETA0(i,j);
            YNY5(i,j)=0;
        end
    else
        ETA5(i,j)=ETA0(i,j);
        YNY5(i,j)=0;
    end
  end
end
% Compute yielding error for markers
if(ynpl>0)
    YERRNOD(iplast)=(ddd/ynpl)^0.5;
end

% Apply criteria for decreasing the timestep
yndt=0;
% Define displacement timestep dtm
dtm=dt;
maxvx=max(max(abs(vx)));
maxvy=max(max(abs(vy)));
if(dtm*maxvx>dxymax*dx)
    dtm=dxymax*dx/maxvx/dtkoefv
    yndt=1;
end
if(dtm*maxvy>dxymax*dy)
    dtm=dxymax*dy/maxvy/dtkoefv
    yndt=1;
end
% Apply maximal stress change condition
dts=dt;
maxDScurrent=max(max(max(abs(DSXX))),max(max(abs(DSXY))));
if(maxDScurrent>DSmax)
    dts=dts/maxDScurrent*DSmax/dtkoefv
    yndt=1;
end
% Apply pressure timestepping condition
dtp=dt;
maxDPcurrent=max(max(abs(DP)));
if(maxDPcurrent>DPmax)
    dtp=dtp/maxDPcurrent*DPmax/dtkoefv
    yndt=1;
end
% Apply thermal timestepping condition
dtt=dt;
maxDTcurrent=max(max(abs(DT)));
if(maxDTcurrent>DTmax)
    dtt=dtt/maxDTcurrent*DTmax/dtkoefv
    yndt=1;
end
% Decrease computational timestep if too many iterations
dtpl=dt;
if(fix((icount+1)/dtstep)*dtstep==(icount+1))
    % Decrease timestep
    dtpl=dtpl/dtkoef
    yndt=1;
end
% Check if further dt reduction possible
if(dt<=dtmin)
    yndt=0;
end

% Reset timestep, restart iteration
ynstop=0;
if(yndt>0)
    dt=max(dtmin,min([dtt, dtm, dts, dtpl, dtp]));
    % Reset old viscoplastic viscosity, temperature
    ETA=ETA00;
    YNY=YNY00;
    tk2=tk1;
elseif(icount>0)
    % Check max P,T change between iterations
    TERRNOD(iplast)=max(max(abs(DT-DT0)));
    PERRNOD(iplast)=max(max(abs(DP-DP0)));
    % Stop iteration
    if(TERRNOD(iplast)<terrmax && PERRNOD(iplast)<perrmax && ...
            (ynpl==0 || iplast==nplast || YERRNOD(iplast)<yerrmax))
        ynstop=1;
    % Repeat iteration
    elseif(ynpl>0)
        % Use new viscoplastic viscosity
        ETA=ETA5;
        YNY=YNY5;
    end
elseif(ynpl>0)
    % Use new viscoplastic viscosity
    ETA=ETA5;
    YNY=YNY5;
end

% Save Last DT and DP
DT0=DT;
DP0=DP;

% Reset/Update iteration counter without dt changes
if(yndt>0)
    icount=0;
else
    icount=icount+1;
end


% Visualise results
if(fix(iplast/visstep)*visstep==iplast || ynstop>0)

figure(3);colormap('Jet');clf
subplot(3,4,1)
pcolor(x/1000,y/1000,log10(ETA0));% caxis([17 21])
shading flat;
axis ij image;
colorbar
title(['log10ETA0, Pa*s timestep=',num2str(timestep)])

subplot(3,4,2)
pcolor(x/1000,y/1000,log10(ETA00));% caxis([17 21])
shading flat;
axis ij image;
colorbar
title(['log10ETA00, Pa*s time=',num2str(timesum),' s'])

subplot(3,4,3)
pcolor(x/1000,y/1000,log10(ETA));% caxis([17 21])
shading flat;
axis ij image;
colorbar
title(['log10ETA, Pa*s iteration=',num2str(iplast)])

subplot(3,4,4)
pcolor(x/1000,y/1000,DSY);% caxis([17 21])
shading flat;
axis ij image;
colorbar
title(['Yield error, Pa for nodes=',num2str(ynpl)])

subplot(3,4,5)
pcolor(xvx/1000,yvx/1000,vx)
shading interp;
axis ij image;
colorbar
title(['Vx, m/s time=',num2str(timesum/(365.25*24*3600)),' Yr'])

subplot(3,4,6)
pcolor(xvy/1000,yvy/1000,vy)
shading interp;
axis ij image;
colorbar
title('Vy, m/s')

subplot(3,4,7)
pcolor(xp/1000,yp/1000,HS)
shading interp;
axis ij image;
colorbar
title('Shear heating, W/m^3')

subplot(3,4,8)
pcolor(xp/1000,yp/1000,log10(EII))
shading interp;
axis ij image;
colorbar
title('logEPSILONii, 1/s')

% Compute/Plot yielding error
subplot(3,1,3)
plot(1:iplast,log10(YERRNOD(1:iplast)),'k');
title(['log Yielding error=',num2str(YERRNOD(iplast)),' Pa for ',num2str(ynpl),' nodes time=',num2str(timesum),' s dt=',num2str(dt),' s'])

pause(.1)
end


% Exiting iteration
if(ynstop>0)
    break
end
% End Plastic iterations on Nodes
end


% Interpolate new viscoplastic viscosity to markers
for m=1:1:marknum
    % Interpolation viscosity from basic nodes
    % Define i,j indexes for the upper left node
    j=fix((xm(m)-x(1))/dx)+1;
    i=fix((ym(m)-y(1))/dy)+1;
    if(j<1)
        j=1;
    elseif(j>Nx-1)
        j=Nx-1;
    end
    if(i<1)
        i=1;
    elseif(i>Ny-1)
        i=Ny-1;
    end
    % Compute distances
    dxmj=xm(m)-x(j);
    dymi=ym(m)-y(i);
    % Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy);
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy);
    wtmi1j1=(dxmj/dx)*(dymi/dy);
    if(YNY(i,j)>0 || YNY(i+1,j)>0 || YNY(i,j+1)>0 || YNY(i+1,j+1)>0)
%         etavpm(m)=ETA(i,j)*wtmij+ETA(i+1,j)*wtmi1j+...
%                 ETA(i,j+1)*wtmij1+ETA(i+1,j+1)*wtmi1j1;
%         etavpm(m)=1/(1/ETA(i,j)*wtmij+1/ETA(i+1,j)*wtmi1j+...
%                 1/ETA(i,j+1)*wtmij1+1/ETA(i+1,j+1)*wtmi1j1);
        etavpm(m)=1/(YNY(i,j)/ETA(i,j)*wtmij+YNY(i+1,j)/ETA(i+1,j)*wtmi1j+...
                YNY(i,j+1)/ETA(i,j+1)*wtmij1+YNY(i+1,j+1)/ETA(i+1,j+1)*wtmi1j1);
        if(etavpm(m)>=etam(tm(m)))
            etavpm(m)=etam(tm(m));
        end
    else
        etavpm(m)=etam(tm(m));
    end
end



% Reset displacement timestep
dtm=dt;


% Apply subgrid stress, pressure diffusion to markers
if(dsubgrids>0)
SXYSUM=zeros(Ny,Nx);
WTSUM=zeros(Ny,Nx);
SXXSUM=zeros(Ny1,Nx1);
PRSUM=zeros(Ny1,Nx1);
WTPSUM=zeros(Ny1,Nx1);
for m=1:1:marknum
    % SIGMA'xx, P
    % Define i,j indexes for the upper left node
    j=fix((xm(m)-xp(1))/dx)+1;
    i=fix((ym(m)-yp(1))/dy)+1;
    if(j<1)
        j=1;
    elseif(j>Nx)
        j=Nx;
    end
    if(i<1)
        i=1;
    elseif(i>Ny)
        i=Ny;
    end
    % Compute distances
    dxmj=xm(m)-xp(j);
    dymi=ym(m)-yp(i);
    % Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy);
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy);
    wtmi1j1=(dxmj/dx)*(dymi/dy);
    % Compute marker-node SIGMA'xx, P difference
    dsxxm0=sxxm(m)-(SXX0(i,j)*wtmij+SXX0(i+1,j)*wtmi1j+...
            SXX0(i,j+1)*wtmij1+SXX0(i+1,j+1)*wtmi1j1);
    dprm0=pr0m(m)-(pr0(i,j)*wtmij+pr0(i+1,j)*wtmi1j+...
            pr0(i,j+1)*wtmij1+pr0(i+1,j+1)*wtmi1j1);
    % Relax stress, P difference
    dsxxm1=dsxxm0*exp(-dsubgrids*dtm/(etam(tm(m))/gggm(tm(m))));
    dprm1=dprm0*exp(-dsubgrids*dtm/(etam(tm(m))/gggm(tm(m))));
    % Correct marker stress, P
    ddsxxm=dsxxm1-dsxxm0;
    sxxm(m)=sxxm(m)+ddsxxm;
    ddprm=dprm1-dprm0;
    pr0m(m)=pr0m(m)+ddprm;
    % Update subgrid diffusion on nodes
    % i,j Node
    SXXSUM(i,j)=SXXSUM(i,j)+ddsxxm*wtmij;
    PRSUM(i,j)=PRSUM(i,j)+ddprm*wtmij;
    WTPSUM(i,j)=WTPSUM(i,j)+wtmij;
    % i+1,j Node
    SXXSUM(i+1,j)=SXXSUM(i+1,j)+ddsxxm*wtmi1j;
    PRSUM(i+1,j)=PRSUM(i+1,j)+ddprm*wtmi1j;
    WTPSUM(i+1,j)=WTPSUM(i+1,j)+wtmi1j;
    % i,j+1 Node
    SXXSUM(i,j+1)=SXXSUM(i,j+1)+ddsxxm*wtmij1;
    PRSUM(i,j+1)=PRSUM(i,j+1)+ddprm*wtmij1;
    WTPSUM(i,j+1)=WTPSUM(i,j+1)+wtmij1;
    % i+1,j+1 Node
    SXXSUM(i+1,j+1)=SXXSUM(i+1,j+1)+ddsxxm*wtmi1j1;
    PRSUM(i+1,j+1)=PRSUM(i+1,j+1)+ddprm*wtmi1j1;
    WTPSUM(i+1,j+1)=WTPSUM(i+1,j+1)+wtmi1j1;

    % SIGMAxy
    % Define i,j indexes for the upper left node
    j=fix((xm(m)-x(1))/dx)+1;
    i=fix((ym(m)-y(1))/dy)+1;
    if(j<1)
        j=1;
    elseif(j>Nx-1)
        j=Nx-1;
    end
    if(i<1)
        i=1;
    elseif(i>Ny-1)
        i=Ny-1;
    end
    % Compute distances
    dxmj=xm(m)-x(j);
    dymi=ym(m)-y(i);
    % Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy);
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy);
    wtmi1j1=(dxmj/dx)*(dymi/dy);
    % Compute marker-node SIGMAxy difference
    dsxym0=sxym(m)-(SXY0(i,j)*wtmij+SXY0(i+1,j)*wtmi1j+...
            SXY0(i,j+1)*wtmij1+SXY0(i+1,j+1)*wtmi1j1);
    % Relax stress difference
    dsxym1=dsxym0*exp(-dsubgrids*dtm/(etam(tm(m))/gggm(tm(m))));
    % Correct marker stress
    ddsxym=dsxym1-dsxym0;
    sxym(m)=sxym(m)+ddsxym;
    % Update subgrid diffusion on nodes
    % i,j Node
    SXYSUM(i,j)=SXYSUM(i,j)+ddsxym*wtmij;
    WTSUM(i,j)=WTSUM(i,j)+wtmij;
    % i+1,j Node
    SXYSUM(i+1,j)=SXYSUM(i+1,j)+ddsxym*wtmi1j;
    WTSUM(i+1,j)=WTSUM(i+1,j)+wtmi1j;
    % i,j+1 Node
    SXYSUM(i,j+1)=SXYSUM(i,j+1)+ddsxym*wtmij1;
    WTSUM(i,j+1)=WTSUM(i,j+1)+wtmij1;
    % i+1,j+1 Node
    SXYSUM(i+1,j+1)=SXYSUM(i+1,j+1)+ddsxym*wtmi1j1;
    WTSUM(i+1,j+1)=WTSUM(i+1,j+1)+wtmi1j1;
end
% Compute DSXXsubgrid, DPsubgrid
DSXXsubgrid=zeros(Ny1,Nx1);
DPsubgrid=zeros(Ny1,Nx1);
% P-nodes
for j=2:1:Nx
    for i=2:1:Ny
        if(WTPSUM(i,j)>0)
            DSXXsubgrid(i,j)=SXXSUM(i,j)/WTPSUM(i,j);
            DPsubgrid(i,j)=PRSUM(i,j)/WTPSUM(i,j);
        end
    end
end
% Correct DSXX, P
DSXX=DSXX-DSXXsubgrid;
DP=DP-DPsubgrid;
% Compute DSXYsubgrid
DSXYsubgrid=zeros(Ny,Nx);
% Basic nodes
for j=1:1:Nx
    for i=1:1:Ny
        if(WTSUM(i,j)>0)
            DSXYsubgrid(i,j)=SXYSUM(i,j)/WTSUM(i,j);
        end
    end
end
% Correct DSXY
DSXY=DSXY-DSXYsubgrid;
end

% Interpolate DSXX,P, DSXY to markers
for m=1:1:marknum
    % SIGMA'xx
    % Define i,j indexes for the upper left node
    j=fix((xm(m)-xp(1))/dx)+1;
    i=fix((ym(m)-yp(1))/dy)+1;
    if(j<2)
        j=2;
    elseif(j>Nx-1)
        j=Nx-1;
    end
    if(i<2)
        i=2;
    elseif(i>Ny-1)
        i=Ny-1;
    end
    % Compute distances
    dxmj=xm(m)-xp(j);
    dymi=ym(m)-yp(i);
    % Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy);
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy);
    wtmi1j1=(dxmj/dx)*(dymi/dy);
    % Update marker by SIGMA'xx, P change 
    sxxm(m)=sxxm(m)+(DSXX(i,j)*wtmij+DSXX(i+1,j)*wtmi1j+...
            DSXX(i,j+1)*wtmij1+DSXX(i+1,j+1)*wtmi1j1);
    pr0m(m)=pr0m(m)+(DP(i,j)*wtmij+DP(i+1,j)*wtmi1j+...
            DP(i,j+1)*wtmij1+DP(i+1,j+1)*wtmi1j1);

    % SIGMAxy
    % Define i,j indexes for the upper left node
    j=fix((xm(m)-x(1))/dx)+1;
    i=fix((ym(m)-y(1))/dy)+1;
    if(j<1)
        j=1;
    elseif(j>Nx-1)
        j=Nx-1;
    end
    if(i<1)
        i=1;
    elseif(i>Ny-1)
        i=Ny-1;
    end
    % Compute distances
    dxmj=xm(m)-x(j);
    dymi=ym(m)-y(i);
    % Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy);
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy);
    wtmi1j1=(dxmj/dx)*(dymi/dy);
    % Update marker by SIGMA'xx change 
    sxym(m)=sxym(m)+(DSXY(i,j)*wtmij+DSXY(i+1,j)*wtmi1j+...
            DSXY(i,j+1)*wtmij1+DSXY(i+1,j+1)*wtmi1j1);
end


% Apply subgrid velocity diffusion on markers
if(dsubgridv>0)
% Vx-nodes
RHOXSUM=zeros(Ny1,Nx1);
VXSUM=zeros(Ny1,Nx1);
% Vy-nodes
RHOYSUM=zeros(Ny1,Nx1);
VYSUM=zeros(Ny1,Nx1);
for m=1:1:marknum
        
    % Interpolation from/to vx-nodes
    % Define i,j indexes for the upper left node
    j=fix((xm(m)-xvx(1))/dx)+1;
    i=fix((ym(m)-yvx(1))/dy)+1;
    if(j<1)
        j=1;
    elseif(j>Nx-1)
        j=Nx-1;
    end
    if(i<1)
        i=1;
    elseif(i>Ny)
        i=Ny;
    end
    % Compute distances
    dxmj=xm(m)-xvx(j);
    dymi=ym(m)-yvx(i);
    % Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy);
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy);
    wtmi1j1=(dxmj/dx)*(dymi/dy);
    % Compute marker-node vx difference
    dvxm0=vx0m(m)-(vx0(i,j)*wtmij+vx0(i+1,j)*wtmi1j+...
            vx0(i,j+1)*wtmij1+vx0(i+1,j+1)*wtmi1j1);
    % Relax vx difference
    dvxm1=dvxm0*exp(-dsubgridv*etavpm(m)/rhom(tm(m))*dtm*(2/dx^2+2/dy^2));
    % Correct marker vx
    ddvxm=dvxm1-dvxm0;
    vx0m(m)=vx0m(m)+ddvxm;
    % Update properties
    % i,j Node
    RHOXSUM(i,j)=RHOXSUM(i,j)+rhom(tm(m))*wtmij;
    VXSUM(i,j)=VXSUM(i,j)+ddvxm*rhom(tm(m))*wtmij;
    % i+1,j Node
    RHOXSUM(i+1,j)=RHOXSUM(i+1,j)+rhom(tm(m))*wtmi1j;
    VXSUM(i+1,j)=VXSUM(i+1,j)+ddvxm*rhom(tm(m))*wtmi1j;
    % i,j+1 Node
    RHOXSUM(i,j+1)=RHOXSUM(i,j+1)+rhom(tm(m))*wtmij1;
    VXSUM(i,j+1)=VXSUM(i,j+1)+ddvxm*rhom(tm(m))*wtmij1;
    % i+1,j+1 Node
    RHOXSUM(i+1,j+1)=RHOXSUM(i+1,j+1)+rhom(tm(m))*wtmi1j1;
    VXSUM(i+1,j+1)=VXSUM(i+1,j+1)+ddvxm*rhom(tm(m))*wtmi1j1;

    % Interpolation from/to vy-nodes
    % Define i,j indexes for the upper left node
    j=fix((xm(m)-xvy(1))/dx)+1;
    i=fix((ym(m)-yvy(1))/dy)+1;
    if(j<1)
        j=1;
    elseif(j>Nx)
        j=Nx;
    end
    if(i<1)
        i=1;
    elseif(i>Ny-1)
        i=Ny-1;
    end
    % Compute distances
    dxmj=xm(m)-xvy(j);
    dymi=ym(m)-yvy(i);
    % Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy);
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy);
    wtmi1j1=(dxmj/dx)*(dymi/dy);
    % Compute marker-node vx difference
    dvym0=vy0m(m)-(vy0(i,j)*wtmij+vy0(i+1,j)*wtmi1j+...
            vy0(i,j+1)*wtmij1+vy0(i+1,j+1)*wtmi1j1);
    % Relax vx difference
    dvym1=dvym0*exp(-dsubgridv*etavpm(m)/rhom(tm(m))*dtm*(2/dx^2+2/dy^2));
    % Correct marker vx
    ddvym=dvym1-dvym0;
    vy0m(m)=vy0m(m)+ddvym;
    % Update properties
    % i,j Node
    RHOYSUM(i,j)=RHOYSUM(i,j)+rhom(tm(m))*wtmij;
    VYSUM(i,j)=VYSUM(i,j)+ddvym*rhom(tm(m))*wtmij;
    % i+1,j Node
    RHOYSUM(i+1,j)=RHOYSUM(i+1,j)+rhom(tm(m))*wtmi1j;
    VYSUM(i+1,j)=VYSUM(i+1,j)+ddvym*rhom(tm(m))*wtmi1j;
    % i,j+1 Node
    RHOYSUM(i,j+1)=RHOYSUM(i,j+1)+rhom(tm(m))*wtmij1;
    VYSUM(i,j+1)=VYSUM(i,j+1)+ddvym*rhom(tm(m))*wtmij1;
    % i+1,j+1 Node
    RHOYSUM(i+1,j+1)=RHOYSUM(i+1,j+1)+rhom(tm(m))*wtmi1j1;
    VYSUM(i+1,j+1)=VYSUM(i+1,j+1)+ddvym*rhom(tm(m))*wtmi1j1;

end
% Compute DVXsubgrid
DVXsubgrid=zeros(Ny1,Nx1);
% Vx-nodes
for j=1:1:Nx1
    for i=1:1:Ny1
        if(RHOXSUM(i,j)>0)
            DVXsubgrid(i,j)=VXSUM(i,j)/RHOXSUM(i,j);
        end
    end
end
% Correct DVX
DVX=DVX-DVXsubgrid;
% Compute DVYsubgrid
DVYsubgrid=zeros(Ny1,Nx1);
% Vy-nodes
for j=1:1:Nx1
    for i=1:1:Ny1
        if(RHOYSUM(i,j)>0)
            DVYsubgrid(i,j)=VYSUM(i,j)/RHOYSUM(i,j);
        end
    end
end
% Correct DVX
DVY=DVY-DVYsubgrid;
end

% Interpolate DVX, DVY to markers
for m=1:1:marknum
    % Interpolation from vx-nodes
    % Define i,j indexes for the upper left node
    j=fix((xm(m)-xvx(1))/dx)+1;
    i=fix((ym(m)-yvx(1))/dy)+1;
    if(j<1)
        j=1;
    elseif(j>Nx-1)
        j=Nx-1;
    end
    if(i<1)
        i=1;
    elseif(i>Ny)
        i=Ny;
    end
    % Compute distances
    dxmj=xm(m)-xvx(j);
    dymi=ym(m)-yvx(i);
    % Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy);
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy);
    wtmi1j1=(dxmj/dx)*(dymi/dy);
    % Update marker vx
    vx0m(m)=vx0m(m)+(DVX(i,j)*wtmij+DVX(i+1,j)*wtmi1j+...
            DVX(i,j+1)*wtmij1+DVX(i+1,j+1)*wtmi1j1);

    % Interpolation from vy-nodes
    % Define i,j indexes for the upper left node
    j=fix((xm(m)-xvy(1))/dx)+1;
    i=fix((ym(m)-yvy(1))/dy)+1;
    if(j<1)
        j=1;
    elseif(j>Nx)
        j=Nx;
    end
    if(i<1)
        i=1;
    elseif(i>Ny-1)
        i=Ny-1;
    end
    % Compute distances
    dxmj=xm(m)-xvy(j);
    dymi=ym(m)-yvy(i);
    % Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy);
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy);
    wtmi1j1=(dxmj/dx)*(dymi/dy);
    % Update marker vy
    vy0m(m)=vy0m(m)+(DVY(i,j)*wtmij+DVY(i+1,j)*wtmi1j+...
            DVY(i,j+1)*wtmij1+DVY(i+1,j+1)*wtmi1j1);
end


% Apply subgrid temperature diffusion on markers
if(dsubgridt>0)
TKSUM=zeros(Ny1,Nx1);
RHOCPSUM=zeros(Ny1,Nx1);
for m=1:1:marknum
    % Define i,j indexes for the upper left node
    j=fix((xm(m)-xp(1))/dx)+1;
    i=fix((ym(m)-yp(1))/dy)+1;
    if(j<1)
        j=1;
    elseif(j>Nx)
        j=Nx;
    end
    if(i<1)
        i=1;
    elseif(i>Ny)
        i=Ny;
    end
    % Compute distances
    dxmj=xm(m)-xp(j);
    dymi=ym(m)-yp(i);
    % Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy);
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy);
    wtmi1j1=(dxmj/dx)*(dymi/dy);
    % Compute marker-node T difference
    dtkm0=tkm(m)-(tk1(i,j)*wtmij+tk1(i+1,j)*wtmi1j+...
            tk1(i,j+1)*wtmij1+tk1(i+1,j+1)*wtmi1j1);
    % Relax temperature difference
    dtkm1=dtkm0*exp(-dsubgridt*km(tm(m))*dtm/rhocpm(tm(m))*(2/dx^2+2/dy^2));
    % Correct marker temperature
    ddtkm=dtkm1-dtkm0;
    tkm(m)=tkm(m)+ddtkm;
    % Update subgrid diffusion on nodes
    % i,j Node
    TKSUM(i,j)=TKSUM(i,j)+ddtkm*rhocpm(tm(m))*wtmij;
    RHOCPSUM(i,j)=RHOCPSUM(i,j)+rhocpm(tm(m))*wtmij;
    % i+1,j Node
    TKSUM(i+1,j)=TKSUM(i+1,j)+ddtkm*rhocpm(tm(m))*wtmi1j;
    RHOCPSUM(i+1,j)=RHOCPSUM(i+1,j)+rhocpm(tm(m))*wtmi1j;
    % i,j+1 Node
    TKSUM(i,j+1)=TKSUM(i,j+1)+ddtkm*rhocpm(tm(m))*wtmij1;
    RHOCPSUM(i,j+1)=RHOCPSUM(i,j+1)+rhocpm(tm(m))*wtmij1;
    % i+1,j+1 Node
    TKSUM(i+1,j+1)=TKSUM(i+1,j+1)+ddtkm*rhocpm(tm(m))*wtmi1j1;
    RHOCPSUM(i+1,j+1)=RHOCPSUM(i+1,j+1)+rhocpm(tm(m))*wtmi1j1;
end
% Compute DTsubgrid
DTsubgrid=zeros(Ny1,Nx1);
% P-nodes
for j=1:1:Nx1
    for i=1:1:Ny1
        if(RHOCPSUM(i,j)>0)
            DTsubgrid(i,j)=TKSUM(i,j)/RHOCPSUM(i,j);
        end
    end
end
% Correct DT
DT=DT-DTsubgrid;
end

% Interpolate DT to markers
for m=1:1:marknum
    % Define i,j indexes for the upper left node
    j=fix((xm(m)-xp(1))/dx)+1;
    i=fix((ym(m)-yp(1))/dy)+1;
    if(j<1)
        j=1;
    elseif(j>Nx)
        j=Nx;
    end
    if(i<1)
        i=1;
    elseif(i>Ny)
        i=Ny;
    end
    % Compute distances
    dxmj=xm(m)-xp(j);
    dymi=ym(m)-yp(i);
    % Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy);
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy);
    wtmi1j1=(dxmj/dx)*(dymi/dy);
    % Update properties
    tkm(m)=tkm(m)+DT(i,j)*wtmij+DT(i+1,j)*wtmi1j+...
            DT(i,j+1)*wtmij1+DT(i+1,j+1)*wtmi1j1;
    % Interpolate tk2 at 1st timestep
    if(timestep==1)
        tkm(m)=tk2(i,j)*wtmij+tk2(i+1,j)*wtmi1j+...
            tk2(i,j+1)*wtmij1+tk2(i+1,j+1)*wtmi1j1;
    end
end

% Compute velocity in pressure nodes
% vx
for j=2:1:Nx
    for i=2:1:Ny
        vxp(i,j)=(vx(i,j)+vx(i,j-1))/2;
    end
end
% Apply BC
% Top
vxp(1,2:Nx-1)=2*vxtop-bctop*vxp(2,2:Nx-1);    
% Bottom
vxp(Ny1,2:Nx-1)=2*vxbottom-bcbottom*vxp(Ny,2:Nx-1);    
% Left
vxp(:,1)=vxp(:,2);
% Right
vxp(:,Nx1)=vxp(:,Nx);
% vy
for j=2:1:Nx
    for i=2:1:Ny
        vyp(i,j)=(vy(i,j)+vy(i-1,j))/2;
    end
end    
% Apply BC
% Left
vyp(2:Ny-1,1)=-bcleft*vyp(2:Ny-1,2);    
% Right
vyp(2:Ny-1,Nx1)=-bcright*vyp(2:Ny-1,Nx);   
% Top
vyp(1,:)=-vyp(2,:);
% Bottom
vyp(Ny1,:)=-vyp(Ny,:);

% Compute rotation rate wyx=1/2(dVy/dx-dVx/dy) for basic nodes
for i=1:1:Ny
    for j=1:1:Nx
        wyx(i,j)=0.5*((vy(i,j+1)-vy(i,j))/dx-(vx(i+1,j)-vx(i,j))/dy);
    end
end


% Compute stress, pressure, velocity state before advection
SXYSUM=zeros(Ny,Nx);
WTSUM=zeros(Ny,Nx);
VXSUM=zeros(Ny1,Nx1);
WTXSUM=zeros(Ny1,Nx1);
VYSUM=zeros(Ny1,Nx1);
WTYSUM=zeros(Ny1,Nx1);
SXXSUM=zeros(Ny1,Nx1);
PRSUM=zeros(Ny1,Nx1);
WTPSUM=zeros(Ny1,Nx1);
for m=1:1:marknum
    % SIGMA'xx, P
    % Define i,j indexes for the upper left node
    j=fix((xm(m)-xp(1))/dx)+1;
    i=fix((ym(m)-yp(1))/dy)+1;
    if(j<1)
        j=1;
    elseif(j>Nx)
        j=Nx;
    end
    if(i<1)
        i=1;
    elseif(i>Ny)
        i=Ny;
    end
    % Compute distances
    dxmj=xm(m)-xp(j);
    dymi=ym(m)-yp(i);
    % Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy);
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy);
    wtmi1j1=(dxmj/dx)*(dymi/dy);
    % Compute nodal SIGMA'xx, P
    sxxm(m)=SXX(i,j)*wtmij+SXX(i+1,j)*wtmi1j+...
            SXX(i,j+1)*wtmij1+SXX(i+1,j+1)*wtmi1j1;
    pr0m(m)=pr(i,j)*wtmij+pr(i+1,j)*wtmi1j+...
            pr(i,j+1)*wtmij1+pr0(i+1,j+1)*wtmi1j1;
    % Update values on nodes
    % i,j Node
    SXXSUM(i,j)=SXXSUM(i,j)+sxxm(m)*wtmij;
    PRSUM(i,j)=PRSUM(i,j)+pr0m(m)*wtmij;
    WTPSUM(i,j)=WTPSUM(i,j)+wtmij;
    % i+1,j Node
    SXXSUM(i+1,j)=SXXSUM(i+1,j)+sxxm(m)*wtmi1j;
    PRSUM(i+1,j)=PRSUM(i+1,j)+pr0m(m)*wtmi1j;
    WTPSUM(i+1,j)=WTPSUM(i+1,j)+wtmi1j;
    % i,j+1 Node
    SXXSUM(i,j+1)=SXXSUM(i,j+1)+sxxm(m)*wtmij1;
    PRSUM(i,j+1)=PRSUM(i,j+1)+pr0m(m)*wtmij1;
    WTPSUM(i,j+1)=WTPSUM(i,j+1)+wtmij1;
    % i+1,j+1 Node
    SXXSUM(i+1,j+1)=SXXSUM(i+1,j+1)+sxxm(m)*wtmi1j1;
    PRSUM(i+1,j+1)=PRSUM(i+1,j+1)+pr0m(m)*wtmi1j1;
    WTPSUM(i+1,j+1)=WTPSUM(i+1,j+1)+wtmi1j1;

    % SIGMAxy
    % Define i,j indexes for the upper left node
    j=fix((xm(m)-x(1))/dx)+1;
    i=fix((ym(m)-y(1))/dy)+1;
    if(j<1)
        j=1;
    elseif(j>Nx-1)
        j=Nx-1;
    end
    if(i<1)
        i=1;
    elseif(i>Ny-1)
        i=Ny-1;
    end
    % Compute distances
    dxmj=xm(m)-x(j);
    dymi=ym(m)-y(i);
    % Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy);
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy);
    wtmi1j1=(dxmj/dx)*(dymi/dy);
    % Compute nodal SIGMAxy
    sxym(m)=SXY(i,j)*wtmij+SXY(i+1,j)*wtmi1j+...
            SXY(i,j+1)*wtmij1+SXY(i+1,j+1)*wtmi1j1;
    % Update value on nodes
    % i,j Node
    SXYSUM(i,j)=SXYSUM(i,j)+sxym(m)*wtmij;
    WTSUM(i,j)=WTSUM(i,j)+wtmij;
    % i+1,j Node
    SXYSUM(i+1,j)=SXYSUM(i+1,j)+sxym(m)*wtmi1j;
    WTSUM(i+1,j)=WTSUM(i+1,j)+wtmi1j;
    % i,j+1 Node
    SXYSUM(i,j+1)=SXYSUM(i,j+1)+sxym(m)*wtmij1;
    WTSUM(i,j+1)=WTSUM(i,j+1)+wtmij1;
    % i+1,j+1 Node
    SXYSUM(i+1,j+1)=SXYSUM(i+1,j+1)+sxym(m)*wtmi1j1;
    WTSUM(i+1,j+1)=WTSUM(i+1,j+1)+wtmi1j1;
    
    % Interpolation from/to vx-nodes
    % Define i,j indexes for the upper left node
    j=fix((xm(m)-xvx(1))/dx)+1;
    i=fix((ym(m)-yvx(1))/dy)+1;
    if(j<1)
        j=1;
    elseif(j>Nx-1)
        j=Nx-1;
    end
    if(i<1)
        i=1;
    elseif(i>Ny)
        i=Ny;
    end
    % Compute distances
    dxmj=xm(m)-xvx(j);
    dymi=ym(m)-yvx(i);
    % Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy);
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy);
    wtmi1j1=(dxmj/dx)*(dymi/dy);
    % Compute nodal vx
    vx0m(m)=vx(i,j)*wtmij+vx(i+1,j)*wtmi1j+...
            vx(i,j+1)*wtmij1+vx(i+1,j+1)*wtmi1j1;
    % Update properties
    % i,j Node
    VXSUM(i,j)=VXSUM(i,j)+vx0m(m)*wtmij;
    WTXSUM(i,j)=WTXSUM(i,j)+wtmij;
    % i+1,j Node
    VXSUM(i+1,j)=VXSUM(i+1,j)+vx0m(m)*wtmi1j;
    WTXSUM(i+1,j)=WTXSUM(i+1,j)+wtmi1j;
    % i,j+1 Node
    VXSUM(i,j+1)=VXSUM(i,j+1)+vx0m(m)*wtmij1;
    WTXSUM(i,j+1)=WTXSUM(i,j+1)+wtmij1;
    % i+1,j+1 Node
    VXSUM(i+1,j+1)=VXSUM(i+1,j+1)+vx0m(m)*wtmi1j1;
    WTXSUM(i+1,j+1)=WTXSUM(i+1,j+1)+wtmi1j1;

    % Interpolation from/to vy-nodes
    % Define i,j indexes for the upper left node
    j=fix((xm(m)-xvy(1))/dx)+1;
    i=fix((ym(m)-yvy(1))/dy)+1;
    if(j<1)
        j=1;
    elseif(j>Nx)
        j=Nx;
    end
    if(i<1)
        i=1;
    elseif(i>Ny-1)
        i=Ny-1;
    end
    % Compute distances
    dxmj=xm(m)-xvy(j);
    dymi=ym(m)-yvy(i);
    % Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy);
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy);
    wtmi1j1=(dxmj/dx)*(dymi/dy);
    % Compute nodal vy
    vy0m(m)=vy(i,j)*wtmij+vy(i+1,j)*wtmi1j+...
            vy(i,j+1)*wtmij1+vy(i+1,j+1)*wtmi1j1;
    % Update properties
    % i,j Node
    VYSUM(i,j)=VYSUM(i,j)+vy0m(m)*wtmij;
    WTYSUM(i,j)=WTYSUM(i,j)+wtmij;
    % i+1,j Node
    VYSUM(i+1,j)=VYSUM(i+1,j)+vy0m(m)*wtmi1j;
    WTYSUM(i+1,j)=WTYSUM(i+1,j)+wtmi1j;
    % i,j+1 Node
    VYSUM(i,j+1)=VYSUM(i,j+1)+vy0m(m)*wtmij1;
    WTYSUM(i,j+1)=WTYSUM(i,j+1)+wtmij1;
    % i+1,j+1 Node
    VYSUM(i+1,j+1)=VYSUM(i+1,j+1)+vy0m(m)*wtmi1j1;
    WTYSUM(i+1,j+1)=WTYSUM(i+1,j+1)+wtmi1j1;
end
% Basic nodes
SXY01=zeros(Ny,Nx);
for j=1:1:Nx
    for i=1:1:Ny
        if(WTSUM(i,j)>0)
            SXY01(i,j)=SXYSUM(i,j)/WTSUM(i,j);
        end
    end
end
% Vx-nodes
VX01=zeros(Ny1,Nx1);
for j=1:1:Nx1
    for i=1:1:Ny1
        if(WTXSUM(i,j)>0)
            VX01(i,j)=VXSUM(i,j)/WTXSUM(i,j);
        end
    end
end
% Vy-nodes
VY01=zeros(Ny1,Nx1);
for j=1:1:Nx1
    for i=1:1:Ny1
        if(WTYSUM(i,j)>0)
            VY01(i,j)=VYSUM(i,j)/WTYSUM(i,j);
        end
    end
end
% P-nodes
SXX01=zeros(Ny1,Nx1);
P01=zeros(Ny1,Nx1);
for j=1:1:Nx1
    for i=1:1:Ny1
        if(WTPSUM(i,j)>0)
            SXX01(i,j)=SXXSUM(i,j)/WTPSUM(i,j);
            P01(i,j)=PRSUM(i,j)/WTPSUM(i,j);
        end
    end
end



% Move markers with 4th order Runge-Kutta
vxm=zeros(4,1);
vym=zeros(4,1);
for m=1:1:marknum
    
    % Interpolate local rotation rate
    % Define i,j indexes for the upper left node
    j=fix((xm(m)-x(1))/dx)+1;
    i=fix((ym(m)-y(1))/dy)+1;
    if(j<1)
        j=1;
    elseif(j>Nx-1)
        j=Nx-1;
    end
    if(i<1)
        i=1;
    elseif(i>Ny-1)
        i=Ny-1;
    end
    % Compute distances
    dxmj=xm(m)-x(j);
    dymi=ym(m)-y(i);
    % Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy);
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy);
    wtmi1j1=(dxmj/dx)*(dymi/dy);
    % Compute vx velocity
    omegam=wyx(i,j)*wtmij+wyx(i+1,j)*wtmi1j+...
    wyx(i,j+1)*wtmij1+wyx(i+1,j+1)*wtmi1j1;
    % Analytical stress rotation using SIGMA'xx=-SIGMA'yy
    THETA=dtm*omegam; % Incremental rotation angle
    sxxmnew=sxxm(m)*cos(THETA)^2-sxxm(m)*sin(THETA)^2-sxym(m)*sin(2*THETA);
    sxymnew=sxxm(m)*sin(2*THETA)+sxym(m)*cos(2*THETA);
    sxxm(m)=sxxmnew; sxym(m)=sxymnew;    
    
    % Save initial marker's coordinates
    xA=xm(m);
    yA=ym(m);
    for rk=1:1:4
        % Interpolate vxp,vyp
        % Define i,j indexes for the upper left node
        j=fix((xm(m)-xp(1))/dx)+1;
        i=fix((ym(m)-yp(1))/dy)+1;
        if(j<1)
            j=1;
        elseif(j>Nx)
            j=Nx;
        end
        if(i<1)
            i=1;
        elseif(i>Ny)
            i=Ny;
        end
        % Compute distances
        dxmj=xm(m)-xp(j);
        dymi=ym(m)-yp(i);
        % Compute weights
        wtmij=(1-dxmj/dx)*(1-dymi/dy);
        wtmi1j=(1-dxmj/dx)*(dymi/dy);    
        wtmij1=(dxmj/dx)*(1-dymi/dy);
        wtmi1j1=(dxmj/dx)*(dymi/dy);
        % Compute vx, vy velocity
        vxm(rk)=vxp(i,j)*wtmij+vxp(i+1,j)*wtmi1j+...
            vxp(i,j+1)*wtmij1+vxp(i+1,j+1)*wtmi1j1;
        vym(rk)=vyp(i,j)*wtmij+vyp(i+1,j)*wtmi1j+...
            vyp(i,j+1)*wtmij1+vyp(i+1,j+1)*wtmi1j1;
        
        % Interpolate vx
        % Define i,j indexes for the upper left node
        j=fix((xm(m)-xvx(1))/dx)+1;
        i=fix((ym(m)-yvx(1))/dy)+1;
        if(j<1)
            j=1;
        elseif(j>Nx-1)
            j=Nx-1;
        end
        if(i<1)
            i=1;
        elseif(i>Ny)
            i=Ny;
        end
        % Compute distances
        dxmj=xm(m)-xvx(j);
        dymi=ym(m)-yvx(i);
        % Compute weights
        wtmij=(1-dxmj/dx)*(1-dymi/dy);
        wtmi1j=(1-dxmj/dx)*(dymi/dy);    
        wtmij1=(dxmj/dx)*(1-dymi/dy);
        wtmi1j1=(dxmj/dx)*(dymi/dy);
        % Compute vx velocity
        vxm(rk)=vpratio*vxm(rk)+(1-vpratio)*(vx(i,j)*wtmij+vx(i+1,j)*wtmi1j+...
            vx(i,j+1)*wtmij1+vx(i+1,j+1)*wtmi1j1);
        
        % Interpolate vy
        % Define i,j indexes for the upper left node
        j=fix((xm(m)-xvy(1))/dx)+1;
        i=fix((ym(m)-yvy(1))/dy)+1;
        if(j<1)
            j=1;
        elseif(j>Nx)
            j=Nx;
        end
        if(i<1)
            i=1;
        elseif(i>Ny-1)
            i=Ny-1;
        end
        % Compute distances
        dxmj=xm(m)-xvy(j);
        dymi=ym(m)-yvy(i);
        % Compute weights
        wtmij=(1-dxmj/dx)*(1-dymi/dy);
        wtmi1j=(1-dxmj/dx)*(dymi/dy);    
        wtmij1=(dxmj/dx)*(1-dymi/dy);
        wtmi1j1=(dxmj/dx)*(dymi/dy);
        % Compute vx velocity
        vym(rk)=vpratio*vym(rk)+(1-vpratio)*(vy(i,j)*wtmij+vy(i+1,j)*wtmi1j+...
            vy(i,j+1)*wtmij1+vy(i+1,j+1)*wtmi1j1);        
        
        % Change coordinates to obtain B,C,D points
        if(rk==1 || rk==2)
            xm(m)=xA+dtm/2*vxm(rk);
            ym(m)=yA+dtm/2*vym(rk);
        elseif(rk==3)
            xm(m)=xA+dtm*vxm(rk);
            ym(m)=yA+dtm*vym(rk);
        end
    end
    % Restore initial coordinates
    xm(m)=xA;
    ym(m)=yA;
    % Compute effective velocity
    vxmeff=1/6*(vxm(1)+2*vxm(2)+2*vxm(3)+vxm(4));
    vymeff=1/6*(vym(1)+2*vym(2)+2*vym(3)+vym(4));
    % Move markers
    xm(m)=xm(m)+dtm*vxmeff;
    ym(m)=ym(m)+dtm*vymeff;
end  


% Compute stress, pressure, velocity state before advection
SXYSUM=zeros(Ny,Nx);
WTSUM=zeros(Ny,Nx);
VXSUM=zeros(Ny1,Nx1);
WTXSUM=zeros(Ny1,Nx1);
VYSUM=zeros(Ny1,Nx1);
WTYSUM=zeros(Ny1,Nx1);
SXXSUM=zeros(Ny1,Nx1);
PRSUM=zeros(Ny1,Nx1);
WTPSUM=zeros(Ny1,Nx1);
for m=1:1:marknum
    % SIGMA'xx, P
    % Define i,j indexes for the upper left node
    j=fix((xm(m)-xp(1))/dx)+1;
    i=fix((ym(m)-yp(1))/dy)+1;
    if(j<1)
        j=1;
    elseif(j>Nx)
        j=Nx;
    end
    if(i<1)
        i=1;
    elseif(i>Ny)
        i=Ny;
    end
    % Compute distances
    dxmj=xm(m)-xp(j);
    dymi=ym(m)-yp(i);
    % Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy);
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy);
    wtmi1j1=(dxmj/dx)*(dymi/dy);
    % Update values on nodes
    % i,j Node
    SXXSUM(i,j)=SXXSUM(i,j)+sxxm(m)*wtmij;
    PRSUM(i,j)=PRSUM(i,j)+pr0m(m)*wtmij;
    WTPSUM(i,j)=WTPSUM(i,j)+wtmij;
    % i+1,j Node
    SXXSUM(i+1,j)=SXXSUM(i+1,j)+sxxm(m)*wtmi1j;
    PRSUM(i+1,j)=PRSUM(i+1,j)+pr0m(m)*wtmi1j;
    WTPSUM(i+1,j)=WTPSUM(i+1,j)+wtmi1j;
    % i,j+1 Node
    SXXSUM(i,j+1)=SXXSUM(i,j+1)+sxxm(m)*wtmij1;
    PRSUM(i,j+1)=PRSUM(i,j+1)+pr0m(m)*wtmij1;
    WTPSUM(i,j+1)=WTPSUM(i,j+1)+wtmij1;
    % i+1,j+1 Node
    SXXSUM(i+1,j+1)=SXXSUM(i+1,j+1)+sxxm(m)*wtmi1j1;
    PRSUM(i+1,j+1)=PRSUM(i+1,j+1)+pr0m(m)*wtmi1j1;
    WTPSUM(i+1,j+1)=WTPSUM(i+1,j+1)+wtmi1j1;

    % SIGMAxy
    % Define i,j indexes for the upper left node
    j=fix((xm(m)-x(1))/dx)+1;
    i=fix((ym(m)-y(1))/dy)+1;
    if(j<1)
        j=1;
    elseif(j>Nx-1)
        j=Nx-1;
    end
    if(i<1)
        i=1;
    elseif(i>Ny-1)
        i=Ny-1;
    end
    % Compute distances
    dxmj=xm(m)-x(j);
    dymi=ym(m)-y(i);
    % Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy);
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy);
    wtmi1j1=(dxmj/dx)*(dymi/dy);
    % Update value on nodes
    % i,j Node
    SXYSUM(i,j)=SXYSUM(i,j)+sxym(m)*wtmij;
    WTSUM(i,j)=WTSUM(i,j)+wtmij;
    % i+1,j Node
    SXYSUM(i+1,j)=SXYSUM(i+1,j)+sxym(m)*wtmi1j;
    WTSUM(i+1,j)=WTSUM(i+1,j)+wtmi1j;
    % i,j+1 Node
    SXYSUM(i,j+1)=SXYSUM(i,j+1)+sxym(m)*wtmij1;
    WTSUM(i,j+1)=WTSUM(i,j+1)+wtmij1;
    % i+1,j+1 Node
    SXYSUM(i+1,j+1)=SXYSUM(i+1,j+1)+sxym(m)*wtmi1j1;
    WTSUM(i+1,j+1)=WTSUM(i+1,j+1)+wtmi1j1;
    
    % Interpolation from vx-nodes
    % Define i,j indexes for the upper left node
    j=fix((xm(m)-xvx(1))/dx)+1;
    i=fix((ym(m)-yvx(1))/dy)+1;
    if(j<1)
        j=1;
    elseif(j>Nx-1)
        j=Nx-1;
    end
    if(i<1)
        i=1;
    elseif(i>Ny)
        i=Ny;
    end
    % Compute distances
    dxmj=xm(m)-xvx(j);
    dymi=ym(m)-yvx(i);
    % Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy);
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy);
    wtmi1j1=(dxmj/dx)*(dymi/dy);
    % Update properties
    % i,j Node
    VXSUM(i,j)=VXSUM(i,j)+vx0m(m)*wtmij;
    WTXSUM(i,j)=WTXSUM(i,j)+wtmij;
    % i+1,j Node
    VXSUM(i+1,j)=VXSUM(i+1,j)+vx0m(m)*wtmi1j;
    WTXSUM(i+1,j)=WTXSUM(i+1,j)+wtmi1j;
    % i,j+1 Node
    VXSUM(i,j+1)=VXSUM(i,j+1)+vx0m(m)*wtmij1;
    WTXSUM(i,j+1)=WTXSUM(i,j+1)+wtmij1;
    % i+1,j+1 Node
    VXSUM(i+1,j+1)=VXSUM(i+1,j+1)+vx0m(m)*wtmi1j1;
    WTXSUM(i+1,j+1)=WTXSUM(i+1,j+1)+wtmi1j1;

    % Interpolation from vy-nodes
    % Define i,j indexes for the upper left node
    j=fix((xm(m)-xvy(1))/dx)+1;
    i=fix((ym(m)-yvy(1))/dy)+1;
    if(j<1)
        j=1;
    elseif(j>Nx)
        j=Nx;
    end
    if(i<1)
        i=1;
    elseif(i>Ny-1)
        i=Ny-1;
    end
    % Compute distances
    dxmj=xm(m)-xvy(j);
    dymi=ym(m)-yvy(i);
    % Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy);
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy);
    wtmi1j1=(dxmj/dx)*(dymi/dy);
    % Update properties
    % i,j Node
    VYSUM(i,j)=VYSUM(i,j)+vy0m(m)*wtmij;
    WTYSUM(i,j)=WTYSUM(i,j)+wtmij;
    % i+1,j Node
    VYSUM(i+1,j)=VYSUM(i+1,j)+vy0m(m)*wtmi1j;
    WTYSUM(i+1,j)=WTYSUM(i+1,j)+wtmi1j;
    % i,j+1 Node
    VYSUM(i,j+1)=VYSUM(i,j+1)+vy0m(m)*wtmij1;
    WTYSUM(i,j+1)=WTYSUM(i,j+1)+wtmij1;
    % i+1,j+1 Node
    VYSUM(i+1,j+1)=VYSUM(i+1,j+1)+vy0m(m)*wtmi1j1;
    WTYSUM(i+1,j+1)=WTYSUM(i+1,j+1)+wtmi1j1;
end
% Basic nodes
SXY02=zeros(Ny,Nx);
for j=1:1:Nx
    for i=1:1:Ny
        if(WTSUM(i,j)>0)
            SXY02(i,j)=SXYSUM(i,j)/WTSUM(i,j);
        end
    end
end
% Vx-nodes
VX02=zeros(Ny1,Nx1);
for j=1:1:Nx1
    for i=1:1:Ny1
        if(WTXSUM(i,j)>0)
            VX02(i,j)=VXSUM(i,j)/WTXSUM(i,j);
        end
    end
end
% Vy-nodes
VY02=zeros(Ny1,Nx1);
for j=1:1:Nx1
    for i=1:1:Ny1
        if(WTYSUM(i,j)>0)
            VY02(i,j)=VYSUM(i,j)/WTYSUM(i,j);
        end
    end
end
% P-nodes
SXX02=zeros(Ny1,Nx1);
P02=zeros(Ny1,Nx1);
for j=1:1:Nx1
    for i=1:1:Ny1
        if(WTPSUM(i,j)>0)
            SXX02(i,j)=SXXSUM(i,j)/WTPSUM(i,j);
            P02(i,j)=PRSUM(i,j)/WTPSUM(i,j);
        end
    end
end



% Correct nodal stresses, pressure, velocity for advection and rotation
SXX0=SXX+SXX02-SXX01;
pr0=pr+P02-P01;
SXY0=SXY+SXY02-SXY01;
vx0=vx+VX02-VX01;
vy0=vy+VY02-VY01;





% Update timesum
timesum=timesum+dtm;
% Trace velocity
velosmax(timestep)=max(max((vxp.^2+vyp.^2).^0.5));
timetrac(timestep)=timesum;
figure(4);colormap('Jet');clf
subplot(2,1,1)
plot(log10(velosmax))
subplot(2,1,2)
plot(timetrac,log10(velosmax))
pause(0.01)

% Visualise results
if(timestep==1 || fix(timestep/visstep)*visstep==timestep)

figure(1);colormap('Jet');clf
subplot(3,4,1)
pcolor(x/1000,y/1000,log10(ETA));% caxis([17 21])
shading flat;
axis ij image;
colorbar
title(['log10ETA, Pa*s timestep=',num2str(timestep)])

subplot(3,4,2)
pcolor(xp/1000,yp/1000,pr)
shading interp;
axis ij image;
colorbar
title('Pressure, Pa')

subplot(3,4,3)
pcolor(xp/1000,yp/1000,vxp)
shading interp;
axis ij image;
colorbar
title(['Vx, m/s time=',num2str(timesum/(365.25*24*3600)),' Yr'])
hold on
quiver(xp(3:5:Nx1)/1000,yp(3:5:Ny1)/1000,vxp(3:5:Ny,3:5:Nx1),vyp(3:5:Ny1,3:5:Nx1),'k')

subplot(3,4,4)
pcolor(xp/1000,yp/1000,vyp)
shading interp;
axis ij image;
colorbar
title('Vy, m/s')
hold on
quiver(xp(3:5:Nx1)/1000,yp(3:5:Ny1)/1000,vxp(3:5:Ny,3:5:Nx1),vyp(3:5:Ny1,3:5:Nx1),'k')

subplot(3,4,5)
pcolor(xp/1000,yp/1000,HS)
shading interp;
axis ij image;
colorbar
title('Shear heating, W/m^3')

subplot(3,4,6)
pcolor(xp/1000,yp/1000,EXX)
shading interp;
axis ij image;
colorbar
title('EPSILONxx, 1/s')

subplot(3,4,7)
pcolor(x/1000,y/1000,EXY)
shading interp;
axis ij image;
colorbar
title('EPSILONxy, 1/s')

subplot(3,4,8)
pcolor(xp/1000,yp/1000,log10(EII))
shading interp;
axis ij image;
colorbar
title('logEPSILONii, 1/s')

subplot(3,4,9)
pcolor(xp/1000,yp/1000,RHO)
shading interp;
axis ij image;
colorbar
title('Density, kg/m^3')

subplot(3,4,10)
pcolor(xp/1000,yp/1000,SXX)
shading interp;
axis ij image;
colorbar
title('SIGMAxx, Pa')

subplot(3,4,11)
pcolor(x/1000,y/1000,SXY)
shading interp;
axis ij image;
colorbar
title('SIGMAxy, Pa')

subplot(3,4,12)
pcolor(xp(2:Nx)/1000,yp(2:Ny)/1000,SII(2:Ny,2:Nx))
shading interp;
axis ij image;
colorbar
title('SIGMAii, Pa')
% Save figure
namefig=['Rate_dependent_' num2str(timestep)];
print ('-djpeg', '-r150',namefig)



% pcolor(x/1000,y/1000,log10(ETA));% caxis([17 21])
% shading flat;
% axis ij image;
% colorbar
% title(['log10ETA, Pa*s timestep=',num2str(timestep),' time=',num2str(timesum/(365.25*24*3600)),' Yr'])


pause(0.1)
end

end
