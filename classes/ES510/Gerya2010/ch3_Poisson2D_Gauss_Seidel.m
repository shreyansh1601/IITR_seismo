% Solution of 2D Poisson equation d2FI/dx2+d2FI/dy2=1 
% with finite differences on a regular grid
% using Gauss-Seidel iteration

% Clean all variables
clear all;
% Clear all figures
clf;

% Numerical model parameters
% Model size, m
xsize=1000; % Horizontal
ysize=1500; % Vertical

% Numbers of nodes
xnum=31; % Horizontal
ynum=41; % Vertical
% Grid step
xstp=xsize/(xnum-1); % Horizontal
ystp=ysize/(ynum-1); % Vertical

% Making vectors for nodal points positions
x=0:xstp:xsize; % Horizontal
y=0:ystp:ysize; % Vertical


% Number of Gauss-Seidel iterations
niter=300;
% Number of iterations before plotting
nplot=10;

% Relaxation factor
TETA=1.5;

% Initial values of PHI
PHI=zeros(ynum,xnum);

% Initial values for right parts of equations
R=zeros(ynum,xnum);
% Define Right Parts
% Process all Grid points
for i=1:1:ynum
  for j=1:1:xnum
    %Boundary nodes
    if(i==1 || i==ynum || j==1 || j==xnum)
        R(i,j)=0;
    else
    %Internal nodes
        R(i,j)=1;
    end
  end
end

% Reset plot iterations counter
ynplot=0;

% Initial values for residuals
dR=zeros(ynum,xnum);

% Iteration Cycle
for k=1:1:niter
    
    % Iterate on internal grid points
    % External ones remain unchanged (FI=0)
    for i=2:1:ynum-1
      for j=2:1:xnum-1
            %d2FI/dx2+d2FI/dy2=1
            %Residual
            dR(i,j)=R(i,j)-((PHI(i,j-1)-2*PHI(i,j)+PHI(i,j+1))/xstp^2+(PHI(i-1,j)-2*PHI(i,j)+PHI(i+1,j))/ystp^2);
            % Compute and assign new value for current unknown
            PHI(i,j)=PHI(i,j)+dR(i,j)/(-2/xstp^2-2/ystp^2)*TETA;
      end
    end

    % Plotting FI and residuals
    % Increase plot iterations counter
    ynplot=ynplot+1;
    if(ynplot==nplot)
        ynplot=0;
        
        % Compute residuals for internal grid points
        for i=2:1:ynum-1
          for j=2:1:xnum-1
                %d2FI/dx2+d2FI/dy2=1
                %Residual
                dR(i,j)=R(i,j)-((PHI(i,j-1)-2*PHI(i,j)+PHI(i,j+1))/xstp^2+(PHI(i-1,j)-2*PHI(i,j)+PHI(i+1,j))/ystp^2);
          end
        end


        % Making new figure
        figure(7);
        %Plotting solution
        subplot(1,2,1);
        surf(x,y,PHI);
        light;
        shading interp;
        colorbar;
        lighting phong;
        %Plotting residuals
        subplot(1,2,2);
        surf(x,y,dR);
        light;
        shading interp;
        colorbar;
        lighting phong;
        pause(0.1);
    end
end
