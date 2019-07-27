% Making manufactured solution benchmark by
% solving Stokes and continuity eq.
% in primitive variable formulation
% with constant viscosity
% using FD with staggered grid
% Analytical solution is given by
% vx=-vx0*sin(2*pi*x/W)*cos(pi*y/H)
% vy=vy0*cos(2*pi*x/W)*sin(pi*y/H)
% P=P0+dP0/dy*y+DP*cos(2*pi*x/W)*sin(pi*y/H)
% Buoyancy therms are then given by x-Stokes and yStokes equations
% RHOgx=-ETA*(d2vx/dx2+d2vx/dy2)+dP/dx
% RHOgy=-ETA*(d2vy/dx2+d2vy/dy2)+dP/dy
% d2vx/dx2=vx0*4*pi^2/W^2*sin(2*pi*x/W)*cos(pi*y/H)
% d2vx/dy2=vx0*pi^2/H^2*sin(2*pi*x/W)*cos(pi*y/H)
% d2vy/dx2=-vy0*4*pi^2/W^2*cos(2*pi*x/W)*sin(pi*y/H)
% d2vy/dy2=-vy0*pi^2/H^2*cos(2*pi*x/W)*sin(pi*y/H)
% dP/dx=-DP*2*pi/W*sin(2*pi*x/W)*sin(pi*y/H)
% dP/dy=dP0/dy+DP*pi/H*cos(2*pi*x/W)*cos(pi*y/H)
% Clearing memory and figures
clear all; clf

% 0) Analytical model parameters
vx0=1e-9;
vy0=3e-9;
P0=1e+5;
dP0dy=0;
DP=1e+6;
H=1500000;
W=1000000;
% 1) Define Numerical model
xsize=W; % Horizontal model size, m
ysize=H; % Vertical model size, m
Nx=81; % Horizontal grid resolution
Ny=101; % Vertical grid resolution
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
pr=zeros(Ny1,Nx1); % Pressure, Pa
vxa=zeros(Ny1,Nx1); % Vx analytical, m/s
vya=zeros(Ny1,Nx1); % Vy analytical, m/s
pra=zeros(Ny1,Nx1); % Pressure analytical, Pa

gy=10; % Gravity acceleration, m/s^2
RHOgx=zeros(Ny1,Nx); % RHO*gx, Pa/m
RHOgy=zeros(Ny,Nx1); % RHO*gy, Pa/m
ETA=1e+19; % Viscosity, Pa*s

% Define Buoyancy terms in internal vx-velocity nodes
for i=2:1:Ny
  for j=2:1:Nx-1
    d2vxdx2=vx0*4*pi^2/W^2*sin(2*pi*xvx(j)/W)*cos(pi*yvx(i)/H);
    d2vxdy2=vx0*pi^2/H^2*sin(2*pi*xvx(j)/W)*cos(pi*yvx(i)/H);
    dPdx=-DP*2*pi/W*sin(2*pi*xvx(j)/W)*sin(pi*yvx(i)/H);
    RHOgx(i,j)=-ETA*(d2vxdx2+d2vxdy2)+dPdx;
  end
end
% Define Buoyancy terms in internal vy-velocity nodes
for i=2:1:Ny-1
  for j=2:1:Nx
    d2vydx2=-vy0*4*pi^2/W^2*cos(2*pi*xvy(j)/W)*sin(pi*yvy(i)/H);
    d2vydy2=-vy0*pi^2/H^2*cos(2*pi*xvy(j)/W)*sin(pi*yvy(i)/H);
    dPdy=dP0dy+DP*pi/H*cos(2*pi*xvy(j)/W)*cos(pi*yvy(i)/H);
    RHOgy(i,j)=-ETA*(d2vydx2+d2vydy2)+dPdy;
  end
end

% Introducing scaled pressure
pscale=ETA/(dx+dy)*2;


% 2) Define global matrixes L(), R()
N=Nx1*Ny1*3; % Global number of unknowns
L=sparse(N,N); % Matrix of coefficients (left part)
R=zeros(N,1); % Vector of right parts

            
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
            % Fictious nodes 1*Vx=0
            if(j==Nx1)
                L(kvx,kvx)=1; % Left part
                R(kvx)=0; % Right part
            end
            % Left boundary
            if(j==1 && i>1 && i<Ny1)
                L(kvx,kvx)=1; % Left part
                R(kvx)=-vx0*sin(2*pi*0/W)*cos(pi*yvx(i)/H); % Right part
            end
            if(j==Nx && i>1 && i<Ny1)
                L(kvx,kvx)=1; % Left part
                R(kvx)=-vx0*sin(2*pi*xsize/W)*cos(pi*yvx(i)/H); % Right part
            end
            
            % Top boundary
            if(i==1 && j<Nx1)
                L(kvx,kvx)=0.5; % Left part
                L(kvx,kvx+3)=0.5; % Left part
                R(kvx)=-vx0*sin(2*pi*xvx(j)/W)*cos(pi*0/H); % Right part
            end
            % Bottom boundary
            if(i==Ny1 && j<Nx1)
                L(kvx,kvx)=0.5; % Left part
                L(kvx,kvx-3)=0.5; % Left part
                R(kvx)=-vx0*sin(2*pi*xvx(j)/W)*cos(pi*ysize/H); % Right part
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
        R(kvx)=-RHOgx(i,j);
        end
        
        % Vy equation External points
        if(j==1 || j==Nx1 || i==1 || i==Ny || i==Ny1)
            % Boundary Condition
            % Fictious nodes
            % 1*Vy=0
            if(i==Ny1)
                L(kvy,kvy)=1; % Left part
                R(kvy)=0; % Right part
            end
            % Top boundary
            if(i==1 && j>1 && j<Nx1)
                L(kvy,kvy)=1; % Left part
                R(kvy)=vy0*cos(2*pi*xvy(j)/W)*sin(pi*0/H); % Right part
            end
            % Bottom boundary
            if(i==Ny && j>1 && j<Nx1)
                L(kvy,kvy)=1; % Left part
                R(kvy)=vy0*cos(2*pi*xvy(j)/W)*sin(pi*ysize/H); % Right part
            end
            % Left boundary
            if(j==1 && i<Ny1)
                L(kvy,kvy)=0.5; % Left part
                L(kvy,kvy+3*Ny1)=0.5; % Left part
                R(kvy)=vy0*cos(2*pi*0/W)*sin(pi*yvy(i)/H); % Right part
            end
            % Right boundary
            if(j==Nx1 && i<Ny1)
                L(kvy,kvy)=0.5; % Left part
                L(kvy,kvy-3*Ny1)=0.5; % Left part
                R(kvy)=vy0*cos(2*pi*xsize/W)*sin(pi*yvy(i)/H); % Right part
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
        R(kvy)=-RHOgy(i,j);
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
                R(kpm)=P0+dP0dy*yp(i)+DP*cos(2*pi*xp(j)/W)*sin(pi*yp(i)/H); % Right part
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
        % Analytical solutions
        % Vx
        if(j<Nx1)
            vxa(i,j)=-vx0*sin(2*pi*xvx(j)/W)*cos(pi*yvx(i)/H);
        end
        % Vy
        if(i<Ny1)
            vya(i,j)=vy0*cos(2*pi*xvy(j)/W)*sin(pi*yvy(i)/H);
        end
        % P
        if(j>1 && j<Nx1 && i>1 && i<Ny1)
            pra(i,j)=P0+dP0dy*yp(i)+DP*cos(2*pi*xp(j)/W)*sin(pi*yp(i)/H);
        end
    end
end


figure(1);colormap('Jet');clf
subplot(3,4,1)
pcolor(xvx,yvx,RHOgx);% caxis([17 21])
shading flat;
axis ij image;
colorbar
title('colormap of RHOgx')

subplot(3,4,5)
pcolor(xvy,yvy,RHOgy);% caxis([17 21])
shading flat;
axis ij image;
colorbar
title('colormap of RHOgy')


subplot(3,4,2)
pcolor(xp(2:Nx),yp(2:Ny),pra(2:Ny,2:Nx))
shading interp;
axis ij image;
colorbar
title('analytical Pressure')

subplot(3,4,6)
pcolor(xvx,yvx,vxa(1:Ny1,1:Nx))
shading interp;
axis ij image;
colorbar
title('analytical of vx')

subplot(3,4,10)
pcolor(xvy,yvy,vya(1:Ny,1:Nx1))
shading interp;
axis ij image;
colorbar
title('analytical vy')


subplot(3,4,3)
pcolor(xp(2:Nx),yp(2:Ny),pr(2:Ny,2:Nx))
shading interp;
axis ij image;
colorbar
title('numerical Pressure')

subplot(3,4,7)
pcolor(xvx,yvx,vx(1:Ny1,1:Nx))
shading interp;
axis ij image;
colorbar
title('numerical vx')

subplot(3,4,11)
pcolor(xvy,yvy,vy(1:Ny,1:Nx1))
shading interp;
axis ij image;
colorbar
title('numerical vy')


subplot(3,4,4)
pcolor(xp(2:Nx),yp(2:Ny),pr(2:Ny,2:Nx)-pra(2:Ny,2:Nx))
shading interp;
axis ij image;
colorbar
title('error Pressure')

subplot(3,4,8)
pcolor(xvx,yvx,vx(1:Ny1,1:Nx)-vxa(1:Ny1,1:Nx))
shading interp;
axis ij image;
colorbar
title('error vx')

subplot(3,4,12)
pcolor(xvy,yvy,vy(1:Ny,1:Nx1)-vya(1:Ny,1:Nx1))
shading interp;
axis ij image;
colorbar
title('error vy')



