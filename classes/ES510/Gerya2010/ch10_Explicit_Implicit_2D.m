% Solution of 2D temperature equation on a regular grid 
% for a non-moving medium with constant conductivity;
% comparison of implicit and explicit method


% Clean all variables
clear all;
% Clear all figures
clf;

% Numerical model parameters
% Boundary conditions = constant temperature
% setup correspond to a squre temperature wave
% Model size, m
xsize=1000000; % Horizontal
ysize=1500000; % Vertical

% Numbers of nodes
xnum=51; % Horizontal
ynum=31; % Vertical
% Grid step
xstp=xsize/(xnum-1); % Horizontal
ystp=ysize/(ynum-1); % Vertical
tnum=20; % number of timesteps

% Material properties
k=3; % Thermal conductivity, W/m/K
cp=1000; % Heat capacity, J/kg/K
rho=3200; % Density, kg/m^3
kappa=k/rho/cp; % Thermal diffusivity, m^2/s

% Making vectors for nodal points positions (basic nodes)
x=0:xstp:xsize; % Horizontal
y=0:ystp:ysize; % Vertical

% Timestep
dtexp=min(xstp,ystp)^2/3/kappa; % Limitation for explicit timestep 
dt=1.0*dtexp; % Timestep, s


% Define initial temperature profile
tback=1000; % background temperature, K
twave=1300; % temperature wave, K
% Creatiing arrays
t0exp=tback*ones(ynum,xnum); % for explicit solving 
t0imp=tback*ones(ynum,xnum); % for implicit solving
for i=1:1:ynum
    for j=1:1:xnum
        % Temperature wave
        if(y(i)>ysize*0.3 && y(i)<ysize*0.7 && x(j)>xsize*0.3 && x(j)<xsize*0.7) 
            t0exp(i,j)=twave; % for explicit solving
            t0imp(i,j)=twave; % for implicit solving
        end
    end
end


% Time cycle
timesum=0; % Elapsed time
for t=1:1:tnum

    % Matrix of coefficients initialization for implicit solving
    L=sparse(xnum*ynum,xnum*ynum);
    % Vector of right part initialization for implicit solving
    R=zeros(xnum*ynum,1);
    
    % Implicit solving of 2D temperature equation:
    % dT/dt=kappa*(d2T/dx2+d2T/dy2)
    % Composing matrix of coefficients L()
    % and vector (column) of right parts R()
    % Process all grid points
    for i=1:1:ynum
        for j=1:1:xnum
            % Global index
            k=(j-1)*ynum+i;
            % Boundary nodes
            if(i==1 || i==ynum || j==1 || j==xnum)
                % Upper boundary
                if(i==1)
                    % Constant temperature: T(i,j)=tback
                    L(k,k)=1;
                    R(k,1)=tback;
%                     % Insulating boundary: T(i,j)=T(i+1,j)
%                     L(k,k)=1;
%                     L(k,k+1)=-1;
%                     R(k,1)=0;                    
                end
                % Lower boundary
                if(i==ynum)
                    % Constant temperature: T(i,j)=tback
                    L(k,k)=1;
                    R(k,1)=tback;
%                     % Insulating boundary: T(i,j)=T(i-1,j)
%                     L(k,k)=1;
%                     L(k,k-1)=-1;
%                     R(k,1)=0;                    
                end
                % Left boundary
                if(j==1 && i>1 && i<ynum)
                    % Constant temperature: T(i,j)=tback
                    L(k,k)=1;
                    R(k,1)=tback;
%                     % Insulating boundary: T(i,j)=T(i,j+1)
%                     L(k,k)=1;
%                     L(k,k+ynum)=-1;
%                     R(k,1)=0;                    
                end
                % Right boundary
                if(j==xnum && i>1 && i<ynum)
                    % Constant temperature: T(i,j)=tback
                    L(k,k)=1;
                    R(k,1)=tback;
%                     % Insulating boundary: T(i,j)=T(i,j-1)
%                     L(k,k)=1;
%                     L(k,k-ynum)=-1;
%                     R(k,1)=0;                    
                end
            % Internal nodes
            else
                % dT/dt=kappa*(d2T/dx2+d2T/dy2)
                % Tnew(i)/dt-kappa*(Tnew(i,j-1)-2*Tnew(i,j)+Tnew(i,j+1))/dx^2-kappa*(Tnew(i-1,j)-2*Tnew(i,j)+Tnew(i+1,j))/dy^2=Told(i,j)/dt
                L(k,k-ynum)=-kappa/xstp^2; % coefficient for T(i,j-1)
                L(k,k+ynum)=-kappa/xstp^2; % coefficient for T(i,j+1)
                L(k,k-1)=-kappa/ystp^2; % coefficient for T(i-1,j)
                L(k,k+1)=-kappa/ystp^2; % coefficient for T(i+1,j)
                L(k,k)=1/dt+2*kappa/xstp^2+2*kappa/ystp^2; % coefficient for T(i,j)
                R(k,1)=t0imp(i,j)/dt;
            end
        end
    end
    % Obtaining solution
    S=L\R;
    % Reloading solution to the new temperature array
    for i=1:1:ynum
        for j=1:1:xnum
            % Global index
            k=(j-1)*ynum+i;
            % Reload solution
            t1imp(i,j)=S(k);
        end
    end
    
    

    % Explicit solving of 2D temperature equation:
    % dT/dt=kappa*(d2T/dx2+d2T/dy2)
    % Process internal grid points
    for i=2:1:ynum-1
        for j=2:1:xnum-1
            % dT/dt=kappa*(d2T/dx2+d2T/dy2)
            % Tnew(i,j)=dt*kappa((Told(i,j-1)-2*Told(i,j)+Told(i,j+1))/dx^2+(Told(i-1,j)-2*Told(i,j)+Told(i+1,j))/dy^2)+Told(i,j)
            t1exp(i,j)=dt*kappa*((t0exp(i,j-1)-2*t0exp(i,j)+t0exp(i,j+1))/xstp^2+(t0exp(i-1,j)-2*t0exp(i,j)+t0exp(i+1,j))/ystp^2)+t0exp(i,j);
        end
    end
    % External points
    % Upper boundary
    % Constant temperature: T(1,j)=tback
    t1exp(1,:)=tback;
%     % Insulating boundary: T(1,j)=T(2,j)
%     t1exp(1,:)=t1exp(2,:);
    % Lower boundary
    % Constant temperature: T(ynum,j)=tback
    t1exp(ynum,:)=tback;
%     % Insulating boundary: T(ynum,j)=T(ynum-1,j)
%     t1exp(ynum,:)=t1exp(ynum-1,:);
    % Left boundary
    % Constant temperature: T(i,1)=tback
    t1exp(2:1:ynum-1,1)=tback;
%     % Insulating boundary: T(i,1)=T(i,2)
%     t1exp(2:1:ynum-1,1)=t1exp(2:1:ynum-1,2);
    % Right boundary
    % Constant temperature: T(i,xnum)=tback
    t1exp(2:1:ynum-1,xnum)=tback;
%     % Insulating boundary: T(i,xnum)=T(i,xnum-1)
%     t1exp(2:1:ynum-1,xnum)=t1exp(2:1:ynum-1,xnum-1);
    
      
    % Open figure
    figure(1);
    % Plotting implicit solution
    subplot(1,2,1);
    pcolor(x/1000,y/1000,t1imp);
    shading interp;
    colorbar;
    axis ij image;
%     caxis([tback-100 twave+100]);
    title(['Implicit: step=',num2str(t),' time,Myr=',num2str(timesum/(1e+6*365.25*24*3600))]);
    % Plotting explicit solution
    subplot(1,2,2);
    pcolor(x/1000,y/1000,t1exp);
    shading interp;
    colorbar;
    axis ij image;
%     caxis([tback-100 twave+100]);
    title(['Explicit: step=',num2str(t),' time,Myr=',num2str(timesum/(1e+6*365.25*24*3600))]);
    
    % Stop for .1 second
    pause(.1);
    
    % Add time counter
    timesum=timesum+dt;
    % Reassign temperature profiles for the next step
    t0exp=t1exp;
    t0imp=t1imp;
end