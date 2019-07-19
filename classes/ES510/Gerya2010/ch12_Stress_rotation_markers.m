% Testing of different stress rotation algoritms
% for 2D stress rotation on advecting markers:
% A. Analytical
% B. Jaumann
% C. Jaumann stress rate with 4th order Runge Kutta
% D. 3D finite angle rotation (Popov et al., 2014)

% Clearing memory and figures
clear all; clf

% 1) Define Numerical model
xsize=1; % Horizontal model size, m
ysize=1; % Vertical model size, m
Nx=51; % Horizontal grid resolution
Ny=51; % Vertical grid resolution
Nx1=Nx+1;
Ny1=Ny+1;
dx=xsize/(Nx-1); % Horizontal grid step, m
dy=ysize/(Ny-1); % Vertical grid step, m
x=0:dx:xsize; % Horizontal coordinates of basic grid points, m
y=0:dy:ysize; % Vertical coordinates of basic grid points, m
xvx=0:dx:xsize; % Horizontal coordinates of vx grid points, m
yvx=-dy/2:dy:ysize+dy/2; % Vertical coordinates of vx grid points, m
xvy=-dx/2:dx:xsize+dx/2; % Horizontal coordinates of vy grid points, m
yvy=0:dy:ysize; % Vertical coordinates of vy grid points, m
xp=-dx/2:dx:xsize+dx/2; % Horizontal coordinates of P grid points, m
yp=-dy/2:dy:ysize+dy/2; % Vertical coordinates of P grid points, m
vx=zeros(Ny1,Nx1); % Vx, m/s
vy=zeros(Ny1,Nx1); % Vy, m/s
vxp=zeros(Ny1,Nx1); % Vx in pressure nodes, m/s
vyp=zeros(Ny1,Nx1); % Vy in pressure nodes, m/s
wyx=zeros(Ny,Nx); % Rotation rate, m/s

% Define markers
Nxm=(Nx-1)*4; % Marker grid resolution in horizontal direction
Nym=(Ny-1)*4; % Marker grid resolution in vertical direction
dxm=xsize/Nxm; % Marker grid step in horizontal direction,m
dym=ysize/Nym; % Marker grid step in vertical direction,m
marknum=Nxm*Nym; % Number of markers
xm=zeros(1,marknum); % Horizontal coordinates, m
ym=zeros(1,marknum); % Vertical coordinates, m
sxxm=zeros(1,marknum); % SIGMAxx', Pa
syym=zeros(1,marknum); % SIGMAyy', Pa
sxym=zeros(1,marknum); % SIGMAxy', Pa

% Compose density array on markers
rp=100000; % Plume radius, m
m=1; % Marker counter
for jm=1:1:Nxm
    for im=1:1:Nym
        % Define marker coordinates
        xm(m)=dxm/2+(jm-1)*dxm+(rand-0.5)*dxm;
        ym(m)=dym/2+(im-1)*dym+(rand-0.5)*dym;
        % Marker properties
        sxxm(m)=1e+6; % SIGMAxx', Pa
        syym(m)=-1e+6; % SIGMAyy', Pa
        sxym(m)=0; % SIGMAxy', Pa
        % Update marker counter
        m=m+1;
    end
end

% Stress rotation parameters
OMEGA=1; % Rotation rate rad/s
LTIME=2*pi/OMEGA; % Time for 1 revolution, s
Nt=100; % Number of timesteps
dt=LTIME/Nt; % Time step size, s
% Defining stress rotation method: 1-4
stressrot=1;
% 1 = A. Analytical
% 2 = B. Jaumann
% 3 = C. Jaumann stress rate with 4th order Runge Kutta
% 4 = D. 3D finite angle rotation (Popov et al., 2014)
% 3D rotation pseudo vector
OMEGAi=[0; 0; 0]; 
% D2. Unit rotation vector 
nu=[0; 0; 0]; 
% 3D stress state
SIGMA3D=[0, 0, 0;
         0, 0, 0;
         0, 0, 0];
% Rotation matrixes 
R1=[1, 0, 0;
    0, 1, 0;
    0, 0, 1];
R2=[0, 0, 0;
    0, 0, 0;
    0, 0, 0];
R3=[0, 0, 0;
    0, 0, 0;
    0, 0, 0];
Rmat=[0, 0, 0;
      0, 0, 0;
      0, 0, 0];

% Select a marker for tracing
mtrace=(((Nx-1)/2)*Ny+2)*4*4+2;
% Arrays to trace stresses and coordinates
SXXT=zeros(1,Nt+1); SXXT(1)=sxxm(mtrace);
SYYT=zeros(1,Nt+1); SYYT(1)=syym(mtrace);
SXYT=zeros(1,Nt+1); SXYT(1)=sxym(mtrace);
XMT=zeros(1,Nt+1); XMT(1)=xm(mtrace);
YMT=zeros(1,Nt+1); YMT(1)=ym(mtrace);

% Define rotational velocity field on the grid
% Vx
for i=1:1:Ny1
    for j=1:1:Nx
        vx(i,j)=OMEGA*(ysize/2-yvx(i));
    end
end
% Vy
for i=1:1:Ny
    for j=1:1:Nx1
        vy(i,j)=OMEGA*(xvy(j)-xsize/2);
    end
end

% Compute rotation rate wyx=1/2(dVy/dx-dVx/dy) for basic nodes
for i=1:1:Ny
    for j=1:1:Nx
        wyx(i,j)=0.5*((vy(i,j+1)-vy(i,j))/dx-(vx(i+1,j)-vx(i,j))/dy);
    end
end

% Compute velocity in pressure nodes
% vxp
for j=2:1:Nx
    for i=1:1:Ny1
        vxp(i,j)=(vx(i,j)+vx(i,j-1))/2;
    end
end
% Apply BC to vxp
% Left
vxp(:,1)=OMEGA*(ysize/2-yp(:));
% Right
vxp(:,Nx1)=OMEGA*(ysize/2-yp(:));
% vy
for j=1:1:Nx1
    for i=2:1:Ny
        vyp(i,j)=(vy(i,j)+vy(i-1,j))/2;
    end
end    
% Apply BC to vyp
% Top
vyp(1,:)=OMEGA*(xp(:)-xsize/2);
% Bottom
vyp(Ny1,:)=OMEGA*(xp(:)-xsize/2);


% Timestepping
vpratio=1/3; % Weight of averaged velocity for moving markers
visstep=10; % Number of steps between visualization
for timestep=1:1:Nt   

% Move markers with 4th order Runge-Kutta
vxm=zeros(4,1);
vym=zeros(4,1);
SXXrate=zeros(4,1);
SYYrate=zeros(4,1);
SXYrate=zeros(4,1);
% Check distance from the model centre
for m=1:1:marknum
if(((xm(m)-xsize/2)^2+(ym(m)-ysize/2)^2)^.5<xsize/2)
    
    % Interpolate local rotation rate wxy
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
    elseif(i>Ny)
        i=Ny;
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
    
    % Rotate stress on markers
    % Update stress state for different methods
    switch(stressrot)
        % A. Analytical
        case 1
        THETA=dt*omegam; % Incremental rotation angle
        sxxmnew=sxxm(m)*cos(THETA)^2+syym(m)*sin(THETA)^2-sxym(m)*sin(2*THETA);
        syymnew=sxxm(m)*sin(THETA)^2+syym(m)*cos(THETA)^2+sxym(m)*sin(2*THETA);
        sxymnew=1/2*(sxxm(m)-syym(m))*sin(2*THETA)+sxym(m)*cos(2*THETA);
        sxxm(m)=sxxmnew; syym(m)=syymnew; sxym(m)=sxymnew;

        % B. Jaumann
        case 2
        THETA=dt*omegam; % Incremental rotation angle
        sxxmnew=sxxm(m)-sxym(m)*2*THETA;
        syymnew=syym(m)+sxym(m)*2*THETA;
        sxymnew=sxym(m)+(sxxm(m)-syym(m))*THETA;
        sxxm(m)=sxxmnew; syym(m)=syymnew; sxym(m)=sxymnew;

        % C. Jaumann stress rate with 4th order Runge Kutta
        case 3
        SXXcur=sxxm(m);
        SYYcur=syym(m);
        SXYcur=sxym(m);
        for rk=1:1:4
            % Compute current Jaumann stress rates
            SXXrate(rk)=-SXYcur*2*omegam;
            SYYrate(rk)=SXYcur*2*omegam;
            SXYrate(rk)=(SXXcur-SYYcur)*omegam;
            % Change current stress state
            if(rk==1 || rk==2)
                SXXcur=sxxm(m)+SXXrate(rk)*dt/2;
                SYYcur=syym(m)+SYYrate(rk)*dt/2;
                SXYcur=sxym(m)+SXYrate(rk)*dt/2;
            elseif(rk==3)
                SXXcur=sxxm(m)+SXXrate(rk)*dt;
                SYYcur=syym(m)+SYYrate(rk)*dt;
                SXYcur=sxym(m)+SXYrate(rk)*dt;
            end
        end
        % Compute effective stress rates
        SXXrateeff=(SXXrate(1)+2*SXXrate(2)+2*SXXrate(3)+SXXrate(4))/6;
        SYYrateeff=(SYYrate(1)+2*SYYrate(2)+2*SYYrate(3)+SYYrate(4))/6;
        SXYrateeff=(SXYrate(1)+2*SXYrate(2)+2*SXYrate(3)+SXYrate(4))/6;
        % Update stresses
        sxxm(m)=sxxm(m)+SXXrateeff*dt;
        syym(m)=syym(m)+SYYrateeff*dt;
        sxym(m)=sxym(m)+SXYrateeff*dt;
        
        % D. 3D finite angle rotation
        case 4
        % 3D rotation pseudo vector
        OMEGAi(3)=omegam; 
        % 3D initial stress state
        SIGMA3D(1,1)=sxxm(m);SIGMA3D(1,2)=sxym(m);
        SIGMA3D(2,1)=sxym(m);SIGMA3D(2,2)=syym(m);
        % D1. Compute vorticity vector magnitude:
        OMEGAmag=(OMEGAi(1)^2+OMEGAi(2)^2+OMEGAi(3)^2)^0.5; 
        % D2. Compute unit rotation vector 
        nu=OMEGAi/OMEGAmag;     
        nx=nu(1);ny=nu(2);nz=nu(3);
        % D3. Integrate incremental rotation angle
        THETA=dt*OMEGAmag;     
        % D4. Evaluate rotation matrix using Euler-Rodrigues formula
        R2(1,2)=-nz;R2(1,3)= ny;
        R2(2,1)= nz;R2(2,3)=-nx;
        R2(3,1)=-ny;R2(3,2)= nx;
        R3(1,1)=nx*nx;R3(1,2)=nx*ny;R3(1,3)=nx*nz;
        R3(2,1)=ny*nx;R3(2,2)=ny*ny;R3(2,3)=ny*nz;
        R3(3,1)=nz*nx;R3(3,2)=nz*ny;R3(3,3)=nz*nz;
        Rmat=cos(THETA)*R1+sin(THETA)*R2+(1-cos(THETA))*R3;
        % D5. Compute rotated stress matrix: SIGMA3D=Rmat*SIGMA3D*Rmat';
        % R3=Rmat*SIGMA3D
        for i=1:1:3
            for j=1:1:3
                R3(i,j)=Rmat(i,1)*SIGMA3D(1,j)+Rmat(i,2)*SIGMA3D(2,j)+Rmat(i,3)*SIGMA3D(3,j);
            end
        end
        % SIGMA3D=R3*Rmat'
        for i=1:1:3
            for j=1:1:3
                SIGMA3D(i,j)=R3(i,1)*Rmat(j,1)+R3(i,2)*Rmat(j,2)+R3(i,3)*Rmat(j,3);
            end
        end
        sxxm(m)=SIGMA3D(1,1);
        syym(m)=SIGMA3D(2,2);
        sxym(m)=SIGMA3D(1,2);
    end
    
    % Advect marker
    % Save initial marker coordinates
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
            xm(m)=xA+dt/2*vxm(rk);
            ym(m)=yA+dt/2*vym(rk);
        elseif(rk==3)
            xm(m)=xA+dt*vxm(rk);
            ym(m)=yA+dt*vym(rk);
        end
    end
    % Restore initial coordinates
    xm(m)=xA;
    ym(m)=yA;
    % Compute effective velocity
    vxmeff=1/6*(vxm(1)+2*vxm(2)+2*vxm(3)+vxm(4));
    vymeff=1/6*(vym(1)+2*vym(2)+2*vym(3)+vym(4));
    % Move markers
    xm(m)=xm(m)+dt*vxmeff;
    ym(m)=ym(m)+dt*vymeff;
end  
end
% Trace stresses and coordinates for selected marker
SXXT(timestep+1)=sxxm(mtrace);
SYYT(timestep+1)=syym(mtrace);
SXYT(timestep+1)=sxym(mtrace);
XMT(timestep+1)=xm(mtrace);
YMT(timestep+1)=ym(mtrace);


% Visualize results
if(fix((timestep)/visstep)*visstep==timestep || timestep==Nt)
figure(1);colormap('Jet');clf
% Vx
subplot(1,2,1)
pcolor(xp,yp,vxp)
shading interp;
axis ij image;
colorbar
title('colormap of vx')
hold on
quiver(xp(3:5:Nx1),yp(3:5:Ny1),vxp(3:5:Ny,3:5:Nx1),vyp(3:5:Ny1,3:5:Nx1),'k')
% Vy
subplot(1,2,2)
pcolor(xp,yp,vyp)
shading interp;
axis ij image;
colorbar
title('colormap of vy')
hold on
quiver(xp(3:5:Nx1),yp(3:5:Ny1),vxp(3:5:Ny,3:5:Nx1),vyp(3:5:Ny1,3:5:Nx1),'k')

figure(2);colormap('Jet');clf
plot(xm,ym,'. k');
hold on
plot(xm(mtrace),ym(mtrace),'o r');
axis ij image


% Plot results
figure(3); clf
% SIGMA'
subplot(2,1,1)
hold on
plot(0:1:timestep,SXXT(1:timestep+1)-SXXT(1),'r');
plot(0:1:timestep,SYYT(1:timestep+1)-SYYT(1),'b');
plot(0:1:timestep,SXYT(1:timestep+1)-SXYT(1),'g');
title('stress changes')
% coordinates
subplot(2,1,2)
hold on
plot(0:1:timestep,XMT(1:timestep+1)-XMT(1),'r');
plot(0:1:timestep,YMT(1:timestep+1)-YMT(1),'b');
title('coordinate changes')

pause(0.1)
end

end

