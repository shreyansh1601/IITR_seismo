% Solving Navier-Stokes and continuity eqs.
% in primitive variable formulation
% with variable viscosity  
% for visco-elasto-plastic RSF material
% using  FD with staggered grid
% and marker in cell techique

% Clearing memory and figures
clear all; clf

% Define Numerical model
xsize=150000; % Horizontal model size, m
ysize=150000; % Vertical model size, m
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
xvx=0:dx:xsize+dx; % Horizontal coordinates of vx grid points, m
yvx=-dy/2:dy:ysize+dy/2; % Vertical coordinates of vx grid points, m
% Vy-nodes
xvy=-dx/2:dx:xsize+dx/2; % Horizontal coordinates of vy grid points, m
yvy=0:dy:ysize+dy; % Vertical coordinates of vy grid points, m
% P-Nodes
xp=-dx/2:dx:xsize+dx/2; % Horizontal coordinates of P grid points, m
yp=-dy/2:dy:ysize+dy/2; % Vertical coordinates of P grid points, m

% Define properties of materials: 
%        Block   Fault   Buffers
rhom    = [2700   2700   2700  ]; % Density, kg/m^3
etam    = [1e+23  1e+23  1e+23 ]; % Viscosity, Pa s
bettam  = [2e-11  2e-11  2e-11 ]; % Compressibility, 1/Pa
gggm    = [3e+10  3e+10  3e+10 ]; % Shear Modulus, Pa
cohesm  = [0e+5   0e+5   0e+5  ]; % Cohesion, Pa
frict0m = [0.2    0.2    0.2   ]; % Reference friction coefficient of RSF
arsfm =   [0.011  0.011  0.011 ]; % a-parameter of RSF
brsfm =   [0.017  0.017  0.001 ]; % b-parameter of RSF
lrsfm =   [0.010  0.010  0.010 ]; % L-parameter of RSF (characteristic slip distance)
omm  =    [20     -1     -1    ]; % Cohesion, Pa
V0=4e-9; % Reference slip velocity of RSF, m/s
DMIN=dy; % Min D-parameter of RSF
DMAX=20000; % Min D-parameter of RSF
KD=1; % Scaling coefficient between THETA*V0 and D=KD*THETA*V0 
% Pressure BC
PCONF=5e+6;



% Nodal arrays
% Basic nodes
ETA=etam(1)*ones(Ny,Nx); % Viscoplastic Viscosity, Pa*s
ETA0=etam(1)*ones(Ny,Nx); % Viscous Viscosity, Pa*s
GGG=gggm(1)*ones(Ny,Nx); % Shear modulus, Pa
COH=cohesm(1)*ones(Ny,Nx); % Cohesion, Pa
FRI0=frict0m(1)*ones(Ny,Nx); % Reference friction coefficient of RSF
OM0=omm(1)*ones(Ny,Nx); % Old state parameter
OM=omm(1)*ones(Ny,Nx); % State parameter
ARSF=arsfm(1)*ones(Ny,Nx); % a-parameter of RSF
BRSF=brsfm(1)*ones(Ny,Nx); % b-parameter of RSF
LRSF=lrsfm(1)*ones(Ny,Nx); % L-parameter of RSF
VPRSF=zeros(Ny,Nx); % Slip velocity of RSF
YNY=zeros(Ny,Nx); % Plastic yielding mark, 1=yes,0=no
EXY=zeros(Ny,Nx); % EPSILONxy, 1/s
SXY=zeros(Ny,Nx); % SIGMAxy, 1/s
SXY0=zeros(Ny,Nx); % SIGMA0xy, 1/s
% Vx-Nodes
RHOX=rhom(1)*ones(Ny1,Nx1); % Density, kg/m^3
vx=zeros(Ny1,Nx1); % vx-velocity m/s
vx0=zeros(Ny1,Nx1); % Old vx-velocity m/s
% Vy-Nodes
RHOY=rhom(1)*ones(Ny1,Nx1); % Density, kg/m^3
vy=zeros(Ny1,Nx1); % vy-velocity m/s
vy0=zeros(Ny1,Nx1); % Old vy-velocity m/s
% P-nodes
BETTA=bettam(1)*ones(Ny1,Nx1); % Compressibility, 1/Pa
ETAP=etam(1)*ones(Ny1,Nx1); % Viscosity, Pa*s
GGGP=gggm(1)*ones(Ny1,Nx1); % Shear modulus, Pa
EXX=zeros(Ny,Nx); % EPSILONxx, 1/s
SXX=zeros(Ny,Nx); % SIGMA'xx, 1/s
SXX0=zeros(Ny,Nx); % SIGMA0'xx, 1/s
vxp=zeros(Ny1,Nx1); % Vx in pressure nodes, m/s
vyp=zeros(Ny1,Nx1); % Vy in pressure nodes, m/s
pr=PCONF*ones(Ny1,Nx1); % Pressure, Pa
pr0=PCONF*ones(Ny1,Nx1); % Old pressure, Pa


% Define Fault
for j=1:1:Nx
    for i=1:1:Ny
        OM0(i,j)=2+log(DMAX/lrsfm(1))-1*exp(-0.5*((y(i)-ysize/2)/10000)^2);
    end
end
OM=OM0;
% DEFF=exp(OM).*LRSF;
% figure(6); clf
% plot(log10(DEFF(:,1)))
% pause(100)

% for j=1:1:Nx
%     if(x(j)>=34000 && x(j)<=116000)
%         % Rate-weakening zone
%         COH(Ny1/2,j)=cohesm(2); % Cohesion, Pa
%         FRI0(Ny1/2,j)=frict0m(2); % Reference friction coefficient of RSF
%         OM0(Ny1/2,j)=omm(2); % Old state parameter
%         OM(Ny1/2,j)=omm(2); % State parameter
%         ARSF(Ny1/2,j)=arsfm(2); % a-parameter of RSF
%         BRSF(Ny1/2,j)=brsfm(2); % b-parameter of RSF
%         LRSF(Ny1/2,j)=lrsfm(2); % L-parameter of RSF
%     elseif(x(j)<=30000 || x(j)>=120000)
%     % Rate-strengthening zone
%         COH(Ny1/2,j)=cohesm(3); % Cohesion, Pa
%         FRI0(Ny1/2,j)=frict0m(3); % Reference friction coefficient of RSF
%         OM0(Ny1/2,j)=omm(3); % Old state parameter
%         OM(Ny1/2,j)=omm(3); % State parameter
%         ARSF(Ny1/2,j)=arsfm(3); % a-parameter of RSF
%         BRSF(Ny1/2,j)=brsfm(3); % b-parameter of RSF
%         LRSF(Ny1/2,j)=lrsfm(3); % L-parameter of RSF
%     elseif(x(j)<34000)
%     % Right transition zone
%         x3=(34000-x(j))/4000;
%         COH(Ny1/2,j)=cohesm(2)+(cohesm(3)-cohesm(2))*x3; % Cohesion, Pa
%         FRI0(Ny1/2,j)=frict0m(2)+(frict0m(3)-frict0m(2))*x3; % Reference friction coefficient of RSF
%         OM0(Ny1/2,j)=omm(2)+(omm(3)-omm(2))*x3; % Old state parameter
%         OM(Ny1/2,j)=omm(2)+(omm(3)-omm(2))*x3; % State parameter
%         ARSF(Ny1/2,j)=arsfm(2)+(arsfm(3)-arsfm(2))*x3; % a-parameter of RSF
%         BRSF(Ny1/2,j)=brsfm(2)+(brsfm(3)-brsfm(2))*x3; % b-parameter of RSF
%         LRSF(Ny1/2,j)=lrsfm(2)+(lrsfm(3)-lrsfm(2))*x3; % L-parameter of RSF
%     elseif(x(j)<120000)
%     % Left transition zone
%         x3=(120000-x(j))/4000;
%         COH(Ny1/2,j)=cohesm(3)+(cohesm(2)-cohesm(3))*x3; % Cohesion, Pa
%         FRI0(Ny1/2,j)=frict0m(3)+(frict0m(2)-frict0m(3))*x3; % Reference friction coefficient of RSF
%         OM0(Ny1/2,j)=omm(3)+(omm(2)-omm(3))*x3; % Old state parameter
%         OM(Ny1/2,j)=omm(3)+(omm(2)-omm(3))*x3; % State parameter
%         ARSF(Ny1/2,j)=arsfm(3)+(arsfm(2)-arsfm(3))*x3; % a-parameter of RSF
%         BRSF(Ny1/2,j)=brsfm(3)+(brsfm(2)-brsfm(3))*x3; % b-parameter of RSF
%         LRSF(Ny1/2,j)=lrsfm(3)+(lrsfm(2)-lrsfm(3))*x3; % L-parameter of RSF
%     end
% end
        
            
% Define global matrixes 
% Mechanical solution: L(), R()
N=Nx1*Ny1*3; % Global number of unknowns
L=sparse(N,N); % Matrix of coefficients (left part)
R=zeros(N,1); % Vector of right parts

% Mechanical boundary conditions: free slip=-1; No Slip=1
bcleft=1;
bcright=1;
bctop=1;
bcbottom=1;
% Top/Bottom velocities
strainrate=1e-13; % Shortening strain rate
vxtop=2e-9;%strainrate*xsize/2;
vxbottom=-2e-9;%-strainrate*xsize/2;

% Thermal boundary conditions: insulation at all boundaries

% Timestepping
dtelastic=1e+8; % Maximal computational timestep, s
dtmin=1e-6; % Minimal computational timestep, s
dt=dtelastic; % Current computational timestep, s
dtkoefv=1.1; % Koefficient to decrese dt for P,T,vx,vy,SIGMA limits 
dtkoef=2; % Koefficient to decrese dt in case of no convergence
dtkoefup=1.1; % Koefficient to increase dt
dtstep=20; % Number of iterations before changing dt
dxymax=0.002; % Max marker movement per time step, grid steps
DSmax=1e+9; % Max stress change per time step, Pa
DPmax=1e+9; % Max pressure change per time step, Pa
DOMmax=0.2; % Max change in RSF state
timesum=0; % Time sum, s
etamin=1e-3; % Lower viscosity cut-off, Pa s
etamax=1e+26; % Upper viscosity cut-off, Pa s
nplast=100000; % Number of plastic iterations
visstep=10; % Periodicity of visualization
savestep=50; % Periodicity of result saving
yerrmax=1e+1; % Tolerance level for yielding error, Pa
perrmax=3e+4; % Tolerance level for pressure error, Pa
YERRNOD=zeros(1,nplast); % Yielding error of nodes
PERRNOD=zeros(1,nplast); % Pressure error of nodes
nsteps=100000; % number of timesteps
namefile='i2elvisSM5_';
timestep=1;
dt00=dt*1.1;
for timestep=timestep:1:nsteps


% Try to increase computational Timestep
if(dt>=dt00)
    dt=min(dt*dtkoefup,dtelastic)
end
dt00=dt;



% Save initial viscoplastic viscosity
ETA00=ETA;
% Save initial yielding nodes
YNY00=YNY;
% Set initial RSF state 
OM=OM0;

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
            R(kpm)=BETTA(i,j)*pr0(i,j)/dt;
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
OM5=OM;
DSY=zeros(Ny,Nx);
ynpl=0;
ddd=0;
dOMdt=zeros(1,4);
dtrobert=dt;
dtrobert1=dt;
dtlapusta=dt;
dtlapusta1=dt;
for i=1:1:Ny
  for j=1:1:Nx
    % Compute second stress invariant
    % SXX, pt are averaged from four surrounding pressure nodes
    SIIB=(SXY(i,j)^2+((SXX(i,j)+SXX(i+1,j)+SXX(i,j+1)+SXX(i+1,j+1))/4)^2)^0.5;
    % Compute pressure
    prB=(pr(i,j)+pr(i+1,j)+pr(i,j+1)+pr(i+1,j+1))/4;
%     prB=PCONF;
        

    % Compute slip velocity for current stress invariant and state
    V=2*V0*sinh(max((SIIB-COH(i,j)),0)/ARSF(i,j)/prB)*...
        exp(-(BRSF(i,j)*OM(i,j)+FRI0(i,j))/ARSF(i,j));
    VPRSF(i,j)=V;
    % Compute new OMEGA
    if(V*dt/LRSF(i,j)>1e-6)
        OM5(i,j)=log(V0/V+(exp(OM0(i,j))-V0/V)*exp(-V*dt/LRSF(i,j)));
    else
        OM5(i,j)=log(exp(OM0(i,j))*(1-V*dt/LRSF(i,j))+V0*dt/LRSF(i,j));
    end
    % Compute yielding stress
    syield=COH(i,j)+prB*ARSF(i,j)*asinh(V/2/V0*...
            exp((BRSF(i,j)*OM5(i,j)+FRI0(i,j))/ARSF(i,j)));
    % Compute new fault thickness
    THETA=exp(OM5(i,j))*LRSF(i,j)/V0;
    D=min(DMAX,max(DMIN,KD*THETA*V0));
    DEFF(i,j)=D;
    % Compute visco-plastic viscosity
    etapl=ETA0(i,j)*syield/(ETA0(i,j)*V/D+syield);
    
    

    % Timestep criterion, Lapusta et al., 2000; Lapusta and Liu, 2009
    B=1/((BETTA(i,j)+BETTA(i+1,j)+BETTA(i,j+1)+BETTA(i+1,j+1))/4);
    vi=(3*B-2*GGG(i,j))/(6*B+2*GGG(i,j));
    k=2/pi*GGG(i,j)/(1-vi)/dx;
    xi=1/4*(k*LRSF(i,j)/ARSF(i,j)/prB-(BRSF(i,j)-...
        ARSF(i,j))/ARSF(i,j))^2-k*LRSF(i,j)/ARSF(i,j)/prB;
    if(xi<0)
        dTETAmax=min(1-(BRSF(i,j)-ARSF(i,j))*prB/(k*LRSF(i,j)),0.2);
    else
        dTETAmax=min(ARSF(i,j)*prB/(k*LRSF(i,j)-(BRSF(i,j)-ARSF(i,j))*prB),0.2);
    end
    dtlapusta=min(dtlapusta,dTETAmax*LRSF(i,j)/V);
    
%     % Herrendoerfer et al. (2018) timestepping:
%     % 1) dt=C/(Max(abs(V/L),abs(V0/L*exp(-phi)))  with C<0.2
%     % 2) t_a=C*etavp/GG with C<0.5 and gg is the shear modulus
%     dtrobert=min(dtrobert,DOMmax/(max(abs(V/LRSF(i,j)),abs(V0/LRSF(i,j)*exp(-OM(i,j))))));
    dtrobert=min(dtrobert,DOMmax*min(LRSF(i,j)/V,LRSF(i,j)/V0*exp(OM(i,j))));
%     dtrobert=min(dtrobert,DOMmax*etapl/GGG(i,j));

    
    
    % Update error for old yielding nodes
    ynn=0;
    if(YNY(i,j)>0)
        ynn=1;
        DSY(i,j)=SIIB-syield;
        ddd=ddd+DSY(i,j)^2;
        ynpl=ynpl+1;
    end
    % Recompute nodal visocity
    ETA5(i,j)=etapl;
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
% Apply state variable timestepping condition
dtom=dt;
maxDOMcurrent=max(max(abs(OM-OM0)));
if(maxDOMcurrent>DOMmax)
    dtom=dtom/maxDOMcurrent*DOMmax/dtkoefv
    yndt=1;
end
if(dtrobert<dt)
    diff=dtrobert-dtlapusta
    yndt=1;
    dtrobert=dtrobert/dtkoefv
end
if(dtlapusta<dt)
    yndt=1;
    dtlapusta=dtlapusta/dtkoefv
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
    dt=max(dtmin,min([dtm, dts, dtpl, dtom, dtrobert, dtlapusta, dtp]));
    % Reset old viscoplastic viscosity, temperature
    ETA=ETA00;
    YNY=YNY00;
    OM=OM0;
elseif(icount>0)
    % Check max P change between iterations
    PERRNOD(iplast)=max(max(abs(DP-DP0)));
    % Stop iteration
    if(PERRNOD(iplast)<perrmax && ...
            (ynpl==0 || iplast==nplast || YERRNOD(iplast)<yerrmax))
        ynstop=1;
    % Repeat iteration
    elseif(ynpl>0)
        % Use new viscoplastic viscosity
        ETA=ETA5;
        YNY=YNY5;
        OM=OM5;
    end
elseif(ynpl>0)
    % Use new viscoplastic viscosity
    ETA=ETA5;
    YNY=YNY5;
    OM=OM5;
end

% Save Last DT and DP
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
% pcolor(x/1000,y/1000,log10(ETA0));% caxis([17 21])
% pcolor(x/1000,y/1000,exp(OM5).*LRSF);% caxis([17 21])
pcolor(x/1000,y/1000,log10(DEFF));% caxis([17 21])
shading flat;
axis ij image;
colorbar
title(['log D, m timestep=',num2str(timestep)])
% title(['log10ETA0, Pa*s timestep=',num2str(timestep)])

subplot(3,4,2)
pcolor(x/1000,y/1000,log10(ETA00));% caxis([17 21])
shading flat;
axis ij image;
colorbar
title(['log10ETA00, Pa*s time=',num2str(timesum/(365.25*24*3600)),' yr'])

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
title(['Vx, m/s time=',num2str(timesum),' s'])

subplot(3,4,6)
pcolor(xvy/1000,yvy/1000,vy)
shading interp;
axis ij image;
colorbar
title('Vy, m/s')

subplot(3,4,7)
pcolor(x/1000,y/1000,OM)
shading interp;
axis ij image;
colorbar
title('State')

subplot(3,4,8)
pcolor(xp/1000,yp/1000,log10(EII))
shading interp;
axis ij image;
colorbar
title('logEPSILONii, 1/s')

% Compute/Plot yielding error
subplot(3,1,3)
plot(1:iplast,log10(YERRNOD(1:iplast)),'k');
title(['log Yielding error=',num2str(YERRNOD(iplast)),' Pa for ',num2str(ynpl),' nodes time=',num2str(timesum/(365.25*24*3600)),' yr dt=',num2str(dt),' s'])

pause(.1)
end


% Exiting iteration
if(ynstop>0)
    break
end
% End Plastic iterations on Nodes
end






% Set old nodal stresses, pressure, state, velocity for advection and rotation
SXX0=SXX;
pr0=pr;
SXY0=SXY;
OM0=OM;
vx0=vx;
vy0=vy;


% Update timesum
timesum=timesum+dtm;
% Trace velocity
% velosmax(timestep)=max(max((vxp.^2+vyp.^2).^0.5));
velosmax(timestep)=max(max(VPRSF));
timetrac(timestep)=timesum;


% Trace slip
if(timestep==1)
    SLIP(timestep,1:Nx)=VPRSF(Ny1/2,1:Nx)*dtm;
else
    SLIP(timestep,1:Nx)=SLIP(timestep-1,1:Nx)+VPRSF(Ny1/2,1:Nx)*dtm;
end


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
pcolor(x/1000,y/1000,OM)
shading interp;
axis ij image;
colorbar
title('State')

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
pcolor(x/1000,y/1000,log10(VPRSF))
shading interp;
axis ij image;
colorbar
title('Vp, m/s')

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
namefig=[namefile num2str(timestep)];
print ('-djpeg', '-r150',namefig)

figure(4);colormap('Jet');clf
subplot(2,1,1)
plot(log10(velosmax))
subplot(2,1,2)
plot(timetrac/(365.25*24*3600),log10(velosmax))
title('Max slip velocity, m/s')

figure(2);colormap('Jet');clf
hold on
plot(x/1000,SLIP(1,:))
for itr=100:100:timestep
plot(x/1000,SLIP(itr,:))
end
title('Slip along the fault, m')

% axis([0 xsize/1000 0 max(SLIP(itr,:))*1.1])

pause(0.1)
end
if(fix(timestep/savestep)*savestep==timestep)
    namemat=[namefile num2str(timestep)];
    save(namemat);
end

end
