% Solution of 2D Stokes and continuity equations with finite differences
% on a regular grid using stream function - vorticity formulation
% for a medium with constant viscosity

% Clean all variables
clear all;
% Clear all figures
clf;

% Numerical model parameters
% Model size, m
xsize=1000000; % Horizontal
ysize=1500000; % Vertical

% Numbers of nodes
xnum=41; % Horizontal
ynum=31; % Vertical
% Grid step
xstp=xsize/(xnum-1); % Horizontal
ystp=ysize/(ynum-1); % Vertical

% Model viscosity
eta=1e+21;

% Gravity acceleration directed downward
gy=10; % m/s^2

% Making vectors for nodal points positions
x=0:xstp:xsize; % Horizontal
y=0:ystp:ysize; % Vertical

% Creating array for density structure (two vertical layers)
rho=zeros(ynum,xnum);
for i=1:1:ynum
  for j=1:1:xnum
    % Horizontal position of the nodal point
    if(x(j)<xsize/2)
        rho(i,j)=3200;  % left layer
    else
        rho(i,j)=3300;  % right layer
    end
  end
end

% Matrix of coefficients initialization
L=sparse(xnum*ynum,xnum*ynum);
% Vector of right part initialization
R=zeros(xnum*ynum,1);

% Solving Poisson equation for worticity
% d2OMEGA/dx2+d2OMEGA/dy2=gy*dRHO/dx
% Composing matrix of coefficients L()
% and vector (column) of right parts R()
% Boundary conditions: OMEGA=0
% Process all Grid points
for i=1:1:ynum
  for j=1:1:xnum
    % Global index for current node
    k=(j-1)*ynum+i;
    %Boundary nodes OMEGA=0
    if(i==1 || i==ynum || j==1 || j==xnum)
        L(k,k)=1;
        R(k,1)=0;
    %Internal nodes
    else
        %Left part: d2OMEGA/dx2+d2OMEGA/dy2
        L(k,k-ynum)=1/xstp^2;
        L(k,k-1)=1/ystp^2;
        L(k,k)=-2/xstp^2-2/ystp^2;
        L(k,k+1)=1/ystp^2;
        L(k,k+ynum)=1/xstp^2;
        % Right part: gy*dRHO/dx
        R(k,1)=gy/eta*(rho(i,j+1)-rho(i,j-1))/2/xstp;
    end
  end
end
   
%Obtaining vector of solutions S()
S=L\R;

% Reload solutions S() to 2D vorticity array OMEGA()
OMEGA=zeros(ynum,xnum);
% Process all Grid points
for i=1:1:ynum
  for j=1:1:xnum
    % Global index for current node
    k=(j-1)*ynum+i;
    OMEGA(i,j)=S(k);
  end
end

% Solving Poisson equation for stream function
% d2PSI/dx2+d2PSI/dy2=OMEGA
% Simplified procedure as below is only possible 
% when boundary conditions for OMEGA and PSI are the same
% othervise i,j cycle is needed as in the previous case
%L=L; % Left parts remain the same
R=S; % Creates right parts from previous solution

%Obtaining vector of solutions S()
S=L\R;

% Reload solutions S() to 2D stream function array PSI()
PSI=zeros(ynum,xnum);
% Process all Grid points
for i=1:1:ynum
  for j=1:1:xnum
    % Global index for current node
    k=(j-1)*ynum+i;
    PSI(i,j)=S(k);
  end
end


% Compute vx,vy for internal nodes
vx=zeros(ynum,xnum);
vy=zeros(ynum,xnum);
% Process internal Grid points
for i=2:1:ynum-1
  for j=2:1:xnum-1
    % vx=dPSI/dy
    vx(i,j)=(PSI(i+1,j)-PSI(i-1,j))/2/xstp;
     % vy=-dPSI/dx
    vy(i,j)=-(PSI(i,j+1)-PSI(i,j-1))/2/xstp;
  end
end


%Plotting solution
% Making new figure
figure(1);

% Plotting vorticity as colormap
subplot(1,2,1);
pcolor(x/1000,y/1000,OMEGA);      % making colormap
shading interp;     % making smooth transitions between colors
colorbar;           % showing colorbar for the map% Plotting Stream function as contours

hold on;            % continuing plotting on the colormap 
[C, h] = contour (x/1000,y/1000,PSI,'k'); % drawing stream funtcion isolines
clabel(C,h,'Color','k');        % labelling the isolines
hold off;           % stop plotting on the colormap

caxis([min(min(OMEGA)) max(max(OMEGA))]); % use color limits for OMEGA
box on;             % making box around the plot
title('worticity, stream function'); % title of the plot
xlabel('x');        % title of horizontal axis
ylabel('y');        % title of vertical axis
axis([0 xsize/1000 0 ysize/1000]); % Making axes limits
axis ij image ;     % direct vertical axis downward, make proper dimensions

% Plotting density as colormap
subplot(1,2,2);
pcolor(x/1000,y/1000,rho);      % making colormap
shading interp;     % making smooth transitions between colors
colorbar;           % showing colorbar for the map

hold on;            % continuing plotting on the colormap 

% Plotting velocity vector as arrows using internal nodes only
quiver(x(2:1:xnum-1)/1000,y(2:1:ynum-1)/1000,vx(2:1:ynum-1,2:1:xnum-1),vy(2:1:ynum-1,2:1:xnum-1),'k'); % making field of arrows

hold off;           % stop plotting on the colormap

box on;             % making box around the plot
title('Density, velocity field');   % title of the plot
xlabel('x');        % title of horizontal axis
ylabel('y');        % title of vertical axis
axis([0 xsize/1000 0 ysize/1000]); % Making axes limits
axis ij image ;     % direct vertical axis downward, make proper dimensions
