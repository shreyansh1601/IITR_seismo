% Solution of 1D temperature equation on a regular grid 
% for a non-moving medium with constant conductivity;
% comparison of implicit and explicit method



% Clean all variables
clear all;
% Clear all figures
clf;

% Boundary conditions = constant temperature
% setup correspond to a squre temperature wave
% Numerical model parameters
xsize=1000; % Model size, m
xnum=51;    % Number of nodes
xstp=xsize/(xnum-1); % Grid step
tnum=50; % number of timesteps
k=3; % Thermal conductivity, W/m/K
cp=1000; % Heat capacity, J/kg/K
rho=3000; % Density, kg/m^3
rhocp=rho*cp; % volumetric heat capacity, j/kg/m^3

kappa=k/rhocp; % Thermal diffusivity, m^2/s

% Timestep
dtexp=xstp^2/2/kappa; % Limitation for explicit timestep
dt=dtexp; % Timestep, s

%Creating vector for nodal point positions
x=0:xstp:xsize;

% Define initial temperature profile
tback=1000; % background temperature, K
twave=1300; % temperature wave, K
for i=1:1:xnum
    % Background profile
    t0exp(i)=tback; % profile for explicit solving
    t0imp(i)=tback; % profile for implicit solving
    % Temperature wave
    if(x(i)>xsize*0.4 && x(i)<xsize*0.6) 
        t0exp(i)=twave; % profile for explicit solving
        t0imp(i)=twave; % profile for implicit solving
    end
end


% Time cycle
timesum=0; % Elapsed time
for t=1:1:tnum

    % Matrix of coefficients initialization for implicit solving
    L=sparse(xnum,xnum);
    % Vector of right part initialization for implicit solving
    R=zeros(xnum,1);
    
    % Implicit solving of 1D temperature equation: dT/dt=kappa*d2T/dx2
    % Composing matrix of coefficients L()
    % and vector (column) of right parts R()
    % First point: T=tback
    L(1,1)=1;
    R(1,1)=tback;
    % Intermediate points
    for i=2:1:xnum-1
        % dT/dt=kappa*d2T/dx2
        % Tnew(i)/dt-kappa*(Tnew(i-1)-2*Tnew(i)+Tnew(i+1))/dx^2=Told(i)/dt
        L(i,i-1)=-kappa/xstp^2;
        L(i,i)=1/dt+2*kappa/xstp^2;
        L(i,i+1)=-kappa/xstp^2;
        R(i,1)=t0imp(i)/dt;
    end
    % Last point:T=tback
    L(xnum,xnum)=1;
    R(xnum,1)=tback;
    % Obtaining solution for implicit temperature profile
    t1imp=L\R;

    % Explicit solving of 1D temperature equation: dT/dt=kappa*d2T/dx2
    % First point: T=tback
    t1exp(1)=tback;
    % Intermediate points
    for i=2:1:xnum-1
        % dT/dt=kappa*d2T/dx2
        % Tnew(i)=dt*kappa(Told(i-1)-2*Told(i)+Told(i+1))/dx^2+Told(i)
        t1exp(i)=dt*kappa*(t0exp(i-1)-2*t0exp(i)+t0exp(i+1))/xstp^2+t0exp(i);
    end
    % Last point:T=tback
    t1exp(xnum)=tback;
    
      
    % Open figure
    figure(1);
    % Plotting implicit solution
    subplot(1,2,1);
    plot(x/1000,t1imp,'k');
    axis([0 xsize/1000 0.9*tback 1.1*twave]);
    title(['Implicit solution, step=',num2str(t)]); 
    % Plotting explicit solution
    subplot(1,2,2);
    plot(x/1000,t1exp,'r');
    axis([0 xsize/1000 0.9*tback 1.1*twave]);
    title(['Explicit solution, step=',num2str(t)]);
    
    % Stop for .1 second
    pause(.1);
    
    % Add time counter
    timesum=timesum+dt;
    % Reassign temperature profiles for the next step
    t0exp=t1exp;
    t0imp=t1imp;
end