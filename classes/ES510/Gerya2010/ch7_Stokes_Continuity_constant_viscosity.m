% Solving Stokes and continuity eq.
% in primitive variable formulation
% with variable viscosity
% using FD with staggered grid

% Clearing memory and figures
clear all; clf

% 1) Define Numerical model
xsize=100000; % Horizontal model size, m
ysize=150000; % Vertical model size, m
Nx=31; % Horizontal grid resolution
Ny=21; % Vertical grid resolution
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
pr=zeros(Ny1,Nx1); % Pressure, Pa
gy=10; % Gravity acceleration, m/s^2
RHO=zeros(Ny,Nx); % Density, kg/m^3
ETA=1e+20; % Viscosity, Pa*s

% Define density in basic nodes
for i=1:1:Ny
  for j=1:1:Nx
    % Horizontal position of the basic nodal point
    if(x(j)<xsize/2)
        % Left layer
        RHO(i,j)=3200;  % density
    else
        % Right layer
        RHO(i,j)=3300;  % density
    end
  end
end

% Introducing scaled pressure
pscale=ETA/(dx+dy)*2;


% 2) Define global matrixes L(), R()
N=Nx1*Ny1*3; % Global number of unknowns
L=sparse(N,N); % Matrix of coefficients (left part)
R=zeros(N,1); % Vector of right parts

% Boundary conditions: free slip=-1; No Slip=1
bcleft=-1;
bcright=-1;
bctop=-1;
bcbottom=-1;

            
% Composing global matrixes L(), R()
% Going through all points of the 2D grid and
% composing respective equations
for j=1:1:Nx1
    for i=1:1:Ny1
        % Define global indexes in algebraic space
        kvx=((j-1)*Ny1+i-1)*3+1; % Vx
        kvy=kvx+1; % Vy
        kpm=kvx+2; % P
        
        % Vx equation External points
        if(i==1 || i==Ny1 || j==1 || j==Nx || j==Nx1)
            % Boundary Condition
            % 1*Vx=0
            L(kvx,kvx)=1; % Left part
            R(kvx)=0; % Right part
            % Top boundary
            if(i==1 && j>1 && j<Nx)
                L(kvx,kvx+3)=bctop; % Left part
            end
            % Bottom boundary
            if(i==Ny1 && j>1 && j<Nx)
                L(kvx,kvx-3)=bcbottom; % Left part
            end
        else
        % Internal points: x-Stokes eq.
        % ETA*(d2Vx/dx^2+d2Vx/dy^2)-dP/dx=0
        %            Vx2
        %             |
        %             |     
        %             |
        %     Vx1-P1-Vx3-P2-Vx5
        %             |
        %             |     
        %             |
        %            Vx4
        %
        % Left part
        L(kvx,kvx-Ny1*3)=ETA/dx^2; % Vx1
        L(kvx,kvx-3)=ETA/dy^2; % Vx2
        L(kvx,kvx)=-2*ETA/dx^2-2*ETA/dy^2; % Vx3
        L(kvx,kvx+3)=ETA/dy^2; % Vx4
        L(kvx,kvx+Ny1*3)=ETA/dx^2; % Vx5
        L(kvx,kpm)=pscale/dx; % P1
        L(kvx,kpm+Ny1*3)=-pscale/dx; % P2
        % Right part
        R(kvx)=0;
        end
        
        % Vy equation External points
        if(j==1 || j==Nx1 || i==1 || i==Ny || i==Ny1)
            % Boundary Condition
            % 1*Vy=0
            L(kvy,kvy)=1; % Left part
            R(kvy)=0; % Right part
            % Left boundary
            if(j==1 && i>1 && i<Ny)
                L(kvy,kvy+3*Ny1)=bcleft; % Left part
            end
            % Right boundary
            if(j==Nx1 && i>1 && i<Ny)
                L(kvy,kvy-3*Ny1)=bcright; % Left part
            end
        else
        % Internal points: y-Stokes eq.
        % ETA*(d2Vy/dx^2+d2Vy/dy^2)-dP/dy=-RHO*gy
        %            Vy2
        %             |
        %             P1    
        %             |
        %     Vy1----Vy3----Vy5
        %             |
        %             P2    
        %             |
        %            Vy4
        %
        % Left part
        L(kvy,kvy-Ny1*3)=ETA/dx^2; % Vy1
        L(kvy,kvy-3)=ETA/dy^2; % Vy2
        L(kvy,kvy)=-2*ETA/dy^2-2*ETA/dx^2; % Vy3
        L(kvy,kvy+3)=ETA/dy^2; % Vy4
        L(kvy,kvy+Ny1*3)=ETA/dx^2; % Vy5
        L(kvy,kpm)=pscale/dy; % P1
        L(kvy,kpm+3)=-pscale/dy; % P2
        
        % Right part
        R(kvy)=-(RHO(i,j-1)+RHO(i,j))/2*gy;
        end
        
        % P equation External points
        if(i==1 || j==1 || i==Ny1 || j==Nx1 ||...
          (i==2 && j==2))
            % Boundary Condition
            % 1*P=0
            L(kpm,kpm)=1; % Left part
            R(kpm)=0; % Right part
            % Real BC
            if(i==2 && j==2)
                L(kpm,kpm)=1*pscale; %Left part
                R(kpm)=0; % Right part
            end
        else
        % Internal points: continuity eq.
        % dVx/dx+dVy/dy=0
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
        % Right part
        R(kpm)=0;
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


figure(1);colormap('Jet');clf
subplot(2,2,1)
pcolor(x,y,RHO);% caxis([17 21])
shading flat;
axis ij image;
colorbar
title('colormap of RHO')
hold on
quiver(xp(3:5:Nx1),yp(3:5:Ny1),vxp(3:5:Ny,3:5:Nx1),vyp(3:5:Ny1,3:5:Nx1),'w')

subplot(2,2,2)
pcolor(xp(2:Nx),yp(2:Ny),pr(2:Ny,2:Nx))
shading interp;
axis ij image;
colorbar
title('colormap of Pressure')

subplot(2,2,3)
pcolor(xvx,yvx,vx(1:Ny1,1:Nx))
shading interp;
axis ij image;
colorbar
title('colormap of vx')

subplot(2,2,4)
pcolor(xvy,yvy,vy(1:Ny,1:Nx1))
shading interp;
axis ij image;
colorbar
title('colormap of vy')


