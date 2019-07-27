% Computation and visualisation of velocity field,
% strain rate, deviatoric strain rate, 
% and second strain rate invariant


% Clearing all variables and arrays
clear all;
% Clearing figures
clf;

% Defining new figure
figure(1)


% Defining model size
xsize=1000000;  % Horizontal size, m
ysize=1500000;  % Vertical size, m
% Defining resolution
xnum=31;    % Horizontal resolution (nodal points)
ynum=31;    % Vertical resolution (nodal points)
% Step between nodal ponts
xstp=xsize/(xnum-1); % Horizontal grid step
ystp=ysize/(ynum-1); % Vertical grid step
% Defining scales for Vx and Vy
vx0=1.1e-9*xsize/2/ysize; % Scale for horizontal velocity, m/s
vy0=1e-9; % Scale for vertical velocity, m/s

% 2D colormap
% Creating a vector of arguments for two axes
x=0:xstp:xsize;
y=0:ystp:ysize;

% Creating velocity function
% Smallest index in an array is 1 (not 0)
for i=1:1:ynum
    for j=1:1:xnum
        % Setup velocity field corresponding 
        % to circulation wit central upwelling
        % Defining Vx = horizontal velocity component distribution
        vx(i,j)=-vx0*sin(pi*x(j)/xsize*2)*cos(pi*y(i)/ysize);
        % Defining Vy = vertical velocity component distribution
        vy(i,j)=vy0*sin(pi*y(i)/ysize)*cos(pi*x(j)/xsize*2);
        
        % Computing partial derivatives 
        % dVx/dx
        dvxdx(i,j)=-vx0*pi/xsize*2*cos(pi*x(j)/xsize*2)*cos(pi*y(i)/ysize);
        % dVx/dy
        dvxdy(i,j)=vx0*pi/ysize*sin(pi*x(j)/xsize*2)*sin(pi*y(i)/ysize);
        % dVy/dy
        dvydy(i,j)=vy0*pi/ysize*cos(pi*y(i)/ysize)*cos(pi*x(j)/xsize*2);
         % dVy/dx
        dvydx(i,j)=-vy0*pi/xsize*2*sin(pi*y(i)/ysize)*sin(pi*x(j)/xsize*2);
        
        % Computing EPSkk=dVx/dx+dVy/dy=div(v)
        epskk(i,j)=dvxdx(i,j)+dvydy(i,j);
        
        % Computing EPS'xx
        eps1xx(i,j)=dvxdx(i,j)-1/3*epskk(i,j);
        % Computing EPS'yy
        eps1yy(i,j)=dvydy(i,j)-1/3*epskk(i,j);
        % Computing EPSxy
        eps1xy(i,j)=1/2*(dvxdy(i,j)+dvydx(i,j));
        % Computing EPSII
        epsii(i,j)=(1/2*(eps1xx(i,j)^2+eps1yy(i,j)^2)+eps1xy(i,j)^2)^0.5;
        
    end
end

figure(1);
% Plotting EPS'xx
subplot(2,3,1);     % defining 2nd plotting area in the 3x1 figure
pcolor(x/1000,y/1000,eps1xx);      % making colormap
box on;             % making box around the plot
title('EPS`xx colormap');   % title of the plot
xlabel('x');        % title of horizontal axis
ylabel('y');        % title of vertical axis
shading interp;     % making smooth transitions between colors
colorbar;           % showing colorbar for the map
axis ij image ;     % direct vertical axis downward, make proper dimensions

% Plotting EPS'yy
subplot(2,3,2);     % defining 2nd plotting area in the 3x1 figure
pcolor(x/1000,y/1000,eps1yy);      % making colormap
box on;             % making box around the plot
title('EPS`yy colormap');   % title of the plot
xlabel('x');        % title of horizontal axis
ylabel('y');        % title of vertical axis
shading interp;     % making smooth transitions between colors
colorbar;           % showing colorbar for the map
axis ij image ;     % direct vertical axis downward, make proper dimensions

% Plotting EPSxy
subplot(2,3,3);     % defining 2nd plotting area in the 3x1 figure
pcolor(x/1000,y/1000,eps1xy);      % making colormap
box on;             % making box around the plot
title('EPSxy colormap');   % title of the plot
xlabel('x');        % title of horizontal axis
ylabel('y');        % title of vertical axis
shading interp;     % making smooth transitions between colors
colorbar;           % showing colorbar for the map
axis ij image ;     % direct vertical axis downward, make proper dimensions

% Plotting velocity vector as arrows
subplot(2,3,4);     % defining 2nd plotting area in the 3x1 figure
quiver(x/1000,y/1000,vx,vy,'k');      % making field of arrows
box on;             % making box around the plot
title('velocity field');   % title of the plot
xlabel('x');        % title of horizontal axis
ylabel('y');        % title of vertical axis
axis([0 xsize/1000 0 ysize/1000]); % Making axes limits
axis ij image ;     % direct vertical axis downward, make proper dimensions


% Plotting EPSII
subplot(2,3,5);     % defining 2nd plotting area in the 3x1 figure
pcolor(x/1000,y/1000,epsii);      % making colormap
box on;             % making box around the plot
title('EPSII colormap');   % title of the plot
xlabel('x');        % title of horizontal axis
ylabel('y');        % title of vertical axis
shading interp;     % making smooth transitions between colors
colorbar;           % showing colorbar for the map
axis ij image ;     % direct vertical axis downward, make proper dimensions

% Plotting EPSkk=div(v)
subplot(2,3,6);     % defining 2nd plotting area in the 3x1 figure
pcolor(x/1000,y/1000,epskk);      % making colormap
box on;             % making box around the plot
title('EPSkk=div(v) colormap');   % title of the plot
xlabel('x');        % title of horizontal axis
ylabel('y');        % title of vertical axis
shading interp;     % making smooth transitions between colors
colorbar;           % showing colorbar for the map
axis ij image ;     % direct vertical axis downward, make proper dimensions
