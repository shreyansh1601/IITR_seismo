% Solution of Poisson equation in 2D with Gauss Seidel iteration
% Density distribution in the model is random

% Clearing all variables and arrays
clear all;
% Clearing figures
clf;

% % Model parameters
% Scaling density
rhoscale=60.0;
% Gravity constant
G=6.67e-11;
% Model size, m
xsize=18000000.0;
ysize=18000000.0;
% Numerical resolution
xnum=100;
ynum=100;
% Grid steps
xstp=xsize./(xnum-1);
ystp=ysize./(ynum-1);
% Poisson equation koefficients
xkf=1/xstp^2;
ykf=1/ystp^2;
xykf=-2/xstp^2-2/ystp^2;

% Iteration Parameters
% Total number of iterations
inum=401;
% Relaxation coefficient for Gauss-Seidel iterations
krelax=0.5;


   
  
% Defining density structure rho()
% Defining initial guesses for gravity potential phi()
% Computing right part of the Poisson equation R()
% Grid points cycle
for i=1:1:ynum;
    for j=1:1:xnum;
        % Gravity potential
        phi(i,j)=0;
        % Define random density
        rho(i,j)=(rand-0.5)*rhoscale;
        % Right part of Poisson equation
        R(i,j)= 10*rho(i,j);
    end
end

% Figures counter
fignum=1;

% Gauss-Seidel iteration circle
for niter=1:1:inum
    
    
    % Plotting Residuals for selected stages
    n=niter-1;
    if(n==0 ||n==5 || n==20 || n==50 || n==100 || n==200)
      % Computing residual for Poisson equation on nodes
        for i=1:1:ynum;
            for j=1:1:xnum;
                % Boundrary conditions 
                if (i==1 | i==ynum | j==1 | j==xnum);
                    % Residuals for boundary conditions
                    residual(i,j)=0;
                    % Solving Poisson equation
                else
                    % Computing residual
                    residual(i,j)=R(i,j)-(xkf*(phi(i,j-1)+phi(i,j+1))+ykf*(phi(i-1,j)+phi(i+1,j))+xykf*phi(i,j));
                end
            end            
        end
        % Plotting Residuals as surface
        figure(1);colormap('Jet');clf
        surf(residual/(10*rhoscale));
%         colormap(gray);
        shading interp;
%         light;
        light('Position',[-1 0 0])
        lighting phong;
        axis tight;
        zlabel('residuals');
%         title(['Solution of Poisson equation, iteration = ',num2str(n)]);
box on
% axis([0 3000 0 100 0 0.010])
ax=gca;
ax.TickDir='out';
ax.LineWidth=1.5;
ax.XTick=[];
ax.YTick=[];
% ax.XTick=[0 1000 2000 3000];
% ax.YTick=[0 50 100];
% ax.ZTick=[0 0.005 0.010];
ax.FontName='Serif';
ax.FontSize=20;
ax.FontWeight='bold';
pause(1)        
print ('-djpeg', '-r300',['smoother',num2str(n)]);
    end
    
    % Solving of the Poisson equation on nodes
    for i=1:1:ynum;
        for j=1:1:xnum;
            % Boundrary conditions 
            if (i==1 | i==ynum | j==1 | j==xnum);
                % Boundary condition for gravity potential
                phi(i,j)=0;
                % Residuals for boundary conditions
                residual(i,j)=0;
            % Solving Poisson equation
            else
                % Computing residual
                resid=R(i,j)-(xkf*(phi(i,j-1)+phi(i,j+1))+ykf*(phi(i-1,j)+phi(i+1,j))+xykf*phi(i,j));
                % Updating solution
                phi(i,j)=phi(i,j)+resid/xykf*krelax;
            end
        end
    end
    
end


