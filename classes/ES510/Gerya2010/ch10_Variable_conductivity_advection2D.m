% Solution of 2D Eulerian temperature equation 
% with advective terms on a regular grid 
% for a moving medium with variable conductivity;
% using of upwind differences for advection of temperature;
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

% Material advection velocity
vx=1e-9; % Horizontal velocity, m/s
vy=1e-9; % Vertical velocity, m/s

% Model properties
% Initial temperature, K 
tback=1000; % background medium,
twave=1300; % hot block
% Thermal conductivity,W/m/K
ktback=3; % background medium 
ktwave=10; % hot block
% Heat capacity, J/K/kg 
cpback=1000; % background medium
cpwave=1100; % hot block
% Density, kg/m^3
rhoback=3200; % background medium
rhowave=3300; % hot block

% Making vectors for nodal points positions (basic nodes)
x=0:xstp:xsize; % Horizontal
y=0:ystp:ysize; % Vertical

% max thermal diffusivity, m^2/s
kappa=max(ktback,ktwave)/min(rhoback,rhowave)/min(cpback,cpwave); 
% Timestep
dtexp=min(xstp,ystp)^2/3/kappa; % Limitation for explicit timestep
% Timestep limitation for advection
if(vx~=0)
    dtexp=min(dtexp,abs(xstp/vx)); % Limitation for horizontal advection timestep
end
if(vy~=0)
    dtexp=min(dtexp,abs(ystp/vy)); % Limitation for vertical advection timestep
end
dt=0.5*dtexp; % Timestep, s

% Creating arrays for Temperature, density, conductivty and heat capacity structure (two vertical layers)
rhocp=rhoback*cpback*ones(ynum,xnum); % volumetric heat capacity, J/K/m^3
kt=ktback*ones(ynum,xnum); % thermal conductivity 
t0exp=tback*ones(ynum,xnum); % for explicit solving 
t0imp=tback*ones(ynum,xnum); % for implicit solving
for i=1:1:ynum
    for j=1:1:xnum
        % Hot block (Temperature wave)
        if(y(i)>ysize*0.3 && y(i)<ysize*0.7 && x(j)>xsize*0.3 && x(j)<xsize*0.7) 
            t0exp(i,j)=twave; % temperature for explicit solving
            t0imp(i,j)=twave; % temperature for implicit solving
            rhocp(i,j)=rhowave*cpwave; % density, kg/m^3
            kt(i,j)=ktwave; % thermal conductivity
        end
    end
end

% Matrix of coefficients initialization for implicit solving
L=sparse(xnum*ynum,xnum*ynum);
% Vector of right part initialization for implicit solving
R=zeros(xnum*ynum,1);


% Time cycle
timesum=0; % Elapsed time
for t=1:1:tnum

    
    % Implicit solving of 2D temperature equation:
    % RHO*Cp*(dT/dt+vx*dT/dx+vy*dT/dy)=d(k*dT/dx)/dx+d(k*dT/dy)/dy
    % Composing matrix of coefficients L()
    % and vector (column) of right parts R()
    % Process all grid pointsd
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
                % RHO*Cp*(dT/dt+vx*dT/dx+vy*dT/dy)=d(k*dT/dx)/dx+d(k*dT/dy)/dy
                % Left part
                % -d(k*dT/dx)/dx
                L(k,k-ynum)=-(kt(i,j-1)+kt(i,j))/2/xstp^2; % coefficient for T(i,j-1)
                L(k,k+ynum)=-(kt(i,j+1)+kt(i,j))/2/xstp^2; % coefficient for T(i,j+1)
                L(k,k)=(kt(i,j-1)+kt(i,j))/2/xstp^2+(kt(i,j+1)+kt(i,j))/2/xstp^2; % coefficient for T(i,j+1)
                % -d(k*dT/dy)/dy
                L(k,k-1)=-(kt(i-1,j)+kt(i,j))/2/ystp^2; % coefficient for T(i-1,j)
                L(k,k+1)=-(kt(i+1,j)+kt(i,j))/2/ystp^2; % coefficient for T(i+1,j)
                L(k,k)=L(k,k)+(kt(i-1,j)+kt(i,j))/2/ystp^2+(kt(i+1,j)+kt(i,j))/2/ystp^2; % ADD coefficient for T(i,j)
                % RHO*Cp*(dT/dt)
                L(k,k)=L(k,k)+rhocp(i,j)/dt; % ADD coefficient for T(i,j)
                % Right part
                R(k,1)=rhocp(i,j)*t0imp(i,j)/dt;
                % Add advective terms to the right part using upwind differences
                % -vx*dT/dx
                if (vx>0)
                    R(k,1)=R(k,1)-rhocp(i,j)*vx*(t0imp(i,j)-t0imp(i,j-1))/xstp;
                else
                    R(k,1)=R(k,1)-rhocp(i,j)*vx*(t0imp(i,j+1)-t0imp(i,j))/xstp;
                end
                % -vy*dT/dy
                if (vy>0)
                    R(k,1)=R(k,1)-rhocp(i,j)*vy*(t0imp(i,j)-t0imp(i-1,j))/ystp;
                else
                    R(k,1)=R(k,1)-rhocp(i,j)*vy*(t0imp(i+1,j)-t0imp(i,j))/ystp;
                end
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
            % RHO*Cp*(dT/dt+vx*dT/dx+vy*dT/dy)=d(k*dT/dx)/dx+d(k*dT/dy)/dy
            % Compute advective terms using upwind differences
            % -vx*dT/dx
            if (vx>0)
                vxdtdx=vx*(t0exp(i,j)-t0exp(i,j-1))/xstp;
            else
                vxdtdx=vx*(t0exp(i,j+1)-t0exp(i,j))/xstp;
            end
            % -vy*dT/dy
            if (vy>0)
                vydtdy=vy*(t0exp(i,j)-t0exp(i-1,j))/ystp;
            else
                vydtdy=vy*(t0exp(i+1,j)-t0exp(i,j))/ystp;
            end
            % Tnew=Told+dt/RHO/Cp*(d(k*dT/dx)/dx+d(k*dT/dy)/dy)
            t1exp(i,j)=t0exp(i,j)+dt/rhocp(i,j)*((kt(i,j+1)+kt(i,j))/2*(t0exp(i,j+1)-t0exp(i,j))-(kt(i,j-1)+kt(i,j))/2*(t0exp(i,j)-t0exp(i,j-1)))/xstp^2;
            t1exp(i,j)=t1exp(i,j)+dt/rhocp(i,j)*((kt(i+1,j)+kt(i,j))/2*(t0exp(i+1,j)-t0exp(i,j))-(kt(i-1,j)+kt(i,j))/2*(t0exp(i,j)-t0exp(i-1,j)))/ystp^2;
            t1exp(i,j)=t1exp(i,j)-dt*(vxdtdx+vydtdy);
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


    % Advect RHO*Cp and k using upwind differences
    % 2D Eulearian advection equation: 
    % dB/dt=-vx*dB/dx-vy*dB/dy
    % Create arrays for new values of RHO*Cp and k
    ktnew=kt;
    rhocpnew=rhocp;
    % Loop through internal nodes
    for i=2:1:ynum-1
        for j=2:1:xnum-1
            % Compute advective terms using upwind differences
            % -vx*dB/dx
            if (vx>0)
                vxdktdx=vx*(kt(i,j)-kt(i,j-1))/xstp;
                vxdrhocpdx=vx*(rhocp(i,j)-rhocp(i,j-1))/xstp;
            else
                vxdktdx=vx*(kt(i,j+1)-kt(i,j))/xstp;
                vxdrhocpdx=vx*(rhocp(i,j+1)-rhocp(i,j))/xstp;
            end
            % -vy*dB/dy
            if (vy>0)
                vydktdy=vy*(kt(i,j)-kt(i-1,j))/ystp;
                vydrhocpdy=vy*(rhocp(i,j)-rhocp(i-1,j))/ystp;
            else
                vydktdy=vy*(kt(i+1,j)-kt(i,j))/ystp;
                vydrhocpdy=vy*(rhocp(i+1,j)-rhocp(i,j))/ystp;
            end
            % Compute new RHO*Cp and k
            rhocpnew(i,j)=rhocp(i,j)-dt*vxdrhocpdx-dt*vydrhocpdy;
            ktnew(i,j)=kt(i,j)-dt*vxdktdx-dt*vydktdy;
        end
    end
    % Update values of RHO*Cp and k
    rhocp=rhocpnew;
    kt=ktnew;
      
    % Open figure
    figure(1);
    % Plotting volumetric heat capacity RHO*Cp
    subplot(2,2,1);
    pcolor(x/1000,y/1000,rhocp);
    shading interp;
    colorbar;
    axis ij image;
    title(['RHO*CP: step=',num2str(t),' time,Myr=',num2str(timesum/(1e+6*365.25*24*3600))]);
    % Plotting thermal conductivity k
    subplot(2,2,2);
    pcolor(x/1000,y/1000,kt);
    shading interp;
    colorbar;
    axis ij image;
    title(['k: step=',num2str(t),' time,Myr=',num2str(timesum/(1e+6*365.25*24*3600))]);
    % Plotting implicit solution
    subplot(2,2,3);
    pcolor(x/1000,y/1000,t1imp);
    shading interp;
    colorbar;
    axis ij image;
%     caxis([tback-100 twave+100]);
    title(['T Implicit: step=',num2str(t),' time,Myr=',num2str(timesum/(1e+6*365.25*24*3600))]);
    % Plotting explicit solution
    subplot(2,2,4);
    pcolor(x/1000,y/1000,t1exp);
    shading interp;
    colorbar;
    axis ij image;
%     caxis([tback-100 twave+100]);
    title(['T Explicit: step=',num2str(t),' time,Myr=',num2str(timesum/(1e+6*365.25*24*3600))]);
    
    % Stop for .1 second
    pause(.1);
    
    % Add time counter
    timesum=timesum+dt;
    % Reassign temperature profiles for the next step
    t0exp=t1exp;
    t0imp=t1imp;
end