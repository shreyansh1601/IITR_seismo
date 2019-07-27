% Testing of different stress rotation algoritms
% for 2D stress rotation:
% A. Analytical
% B. Jaumann
% C. Jaumann stress rate with 4th order Runge Kutta
% D. 3D finite angle rotation (Popov et al., 2014)


% Clearing memory and figures
clear all; clf

% Define numerical model
OMEGA=1; % Rotation rate rad/s
LTIME=2*pi/OMEGA; % Time for 1 revolution, s
Nt=100; % Number of timesteps
dt=LTIME/Nt; % Time step size, s
SXX0=1e+6; % Initial SIGMA'xx stress component, Pa
SYY0=-1e+6; % Initial SIGMA'yy stress component, Pa
SXY0=0; % Initial SIGMA'xx stress component, Pa
% 3D rotation pseudo vector
OMEGAi=[0;     0;    OMEGA]; 
% 3D initial stress state
SIGMA0=[SXX0,  SXY0, 0;
        SXY0,  SYY0, 0;
        0,     0,    0];
% Set initial stress state for different methods
% A. Analytical
SXXA=zeros(1,Nt+1); SXXA(1)=SXX0;
SYYA=zeros(1,Nt+1); SYYA(1)=SYY0;
SXYA=zeros(1,Nt+1); SXYA(1)=SXY0;
% B. Jaumann
SXXJ=zeros(1,Nt+1); SXXJ(1)=SXX0;
SYYJ=zeros(1,Nt+1); SYYJ(1)=SYY0;
SXYJ=zeros(1,Nt+1); SXYJ(1)=SXY0;
% C. Jaumann stress rate with 4th order Runge Kutta
SXXJRK=zeros(1,Nt+1); SXXJRK(1)=SXX0;
SYYJRK=zeros(1,Nt+1); SYYJRK(1)=SYY0;
SXYJRK=zeros(1,Nt+1); SXYJRK(1)=SXY0;
SXXrate=zeros(4,1);
SYYrate=zeros(4,1);
SXYrate=zeros(4,1);
% C. 3D finite angle rotation
SIGMA3D=SIGMA0;
SXX3D=zeros(1,Nt+1); SXX3D(1)=SXX0;
SYY3D=zeros(1,Nt+1); SYY3D(1)=SYY0;
SXY3D=zeros(1,Nt+1); SXY3D(1)=SXY0;

% Timestepping
for t=1:1:Nt
    % Update stress state for different methods
    % A. Analytical
    THETA=dt*OMEGA; % Incremental rotation angle
    SXXA(t+1)=SXXA(t)*cos(THETA)^2+SYYA(t)*sin(THETA)^2-SXYA(t)*sin(2*THETA);
    SYYA(t+1)=SXXA(t)*sin(THETA)^2+SYYA(t)*cos(THETA)^2+SXYA(t)*sin(2*THETA);
    SXYA(t+1)=1/2*(SXXA(t)-SYYA(t))*sin(2*THETA)+SXYA(t)*cos(2*THETA);
    
    % B. Jaumann
    THETA=dt*OMEGA; % Incremental rotation angle
    SXXJ(t+1)=SXXJ(t)-SXYJ(t)*2*THETA;
    SYYJ(t+1)=SYYJ(t)+SXYJ(t)*2*THETA;
    SXYJ(t+1)=SXYJ(t)+(SXXJ(t)-SYYJ(t))*THETA;
    
    % C. Jaumann stress rate with 4th order Runge Kutta
    SXXcur=SXXJRK(t);
    SYYcur=SYYJRK(t);
    SXYcur=SXYJRK(t);
    for rk=1:1:4
        % Compute current Jaumann stress rates
        SXXrate(rk)=-SXYcur*2*OMEGA;
        SYYrate(rk)=SXYcur*2*OMEGA;
        SXYrate(rk)=(SXXcur-SYYcur)*OMEGA;
        % Change current stress state
        if(rk==1 || rk==2)
            SXXcur=SXXJRK(t)+SXXrate(rk)*dt/2;
            SYYcur=SYYJRK(t)+SYYrate(rk)*dt/2;
            SXYcur=SXYJRK(t)+SXYrate(rk)*dt/2;
        elseif(rk==3)
            SXXcur=SXXJRK(t)+SXXrate(rk)*dt;
            SYYcur=SYYJRK(t)+SYYrate(rk)*dt;
            SXYcur=SXYJRK(t)+SXYrate(rk)*dt;
        end
    end
    % Compute effective stress rates
    SXXrateeff=(SXXrate(1)+2*SXXrate(2)+2*SXXrate(3)+SXXrate(4))/6;
    SYYrateeff=(SYYrate(1)+2*SYYrate(2)+2*SYYrate(3)+SYYrate(4))/6;
    SXYrateeff=(SXYrate(1)+2*SXYrate(2)+2*SXYrate(3)+SXYrate(4))/6;
    % Update stresses
    SXXJRK(t+1)=SXXJRK(t)+SXXrateeff*dt;
    SYYJRK(t+1)=SYYJRK(t)+SYYrateeff*dt;
    SXYJRK(t+1)=SXYJRK(t)+SXYrateeff*dt;
    
    % D. 3D finite angle rotation
    % D1. Compute vorticity vector magnitude:
    OMEGAmag=norm(OMEGAi); 
    % D2. Compute unit rotation vector 
    n=OMEGAi/OMEGAmag;     
    nx=n(1); 
    ny=n(2); 
    nz=n(3);
    % D3. Integrate incremental rotation angle
    THETA = dt*OMEGAmag;     
    % D4. Evaluate rotation matrix using Euler-Rodrigues formula
    R1=[ 1,    0,     0;
         0,    1,     0;
         0,    0,     1];
    R2=[ 0,   -nz,    ny;
         nz,   0,    -nx;
        -ny,   nx,    0 ];
    R3=[nx*nx, nx*ny, nx*nz;
        ny*nx, ny*ny, ny*nz;
        nz*nx, nz*ny, nz*nz];
    Rmat=cos(THETA)*R1+sin(THETA)*R2+(1-cos(THETA))*R3;
    % D5. Compute rotated stress matrix 
    SIGMA3D=Rmat*SIGMA3D*Rmat';
    SXX3D(t+1)=SIGMA3D(1,1);
    SYY3D(t+1)=SIGMA3D(2,2);
    SXY3D(t+1)=SIGMA3D(1,2);
end

% Plot results
figure(1); clf
% SIGMA'xx
subplot(3,1,1)
hold on
plot(0:1:Nt,SXXA,'r');
plot(0:1:Nt,SXXJ,'b');
plot(0:1:Nt,SXXJRK,'g');
plot(0:1:Nt,SXX3D,'.k');
title('SIGMAxx')
% SIGMAyy'
subplot(3,1,2)
hold on
plot(0:1:Nt,SYYA,'r');
plot(0:1:Nt,SYYJ,'b');
plot(0:1:Nt,SYYJRK,'g');
plot(0:1:Nt,SYY3D,'.k');
title('SIGMAyy')
% SIGMAxy'
subplot(3,1,3)
hold on
plot(0:1:Nt,SXYA,'r');
plot(0:1:Nt,SXYJ,'b');
plot(0:1:Nt,SXYJRK,'g');
plot(0:1:Nt,SXY3D,'.k');
title('SIGMAxy')

% Plot error
figure(2); clf
% SIGMA'xx
subplot(3,1,1)
hold on
plot(0:1:Nt,SXXJ-SXXA,'.b');
% plot(0:1:Nt,SXXJRK-SXXA,'g');
% plot(0:1:Nt,SXX3D-SXXA,'.k');
title('SIGMAxx')
% SIGMAyy'
subplot(3,1,2)
hold on
% plot(0:1:Nt,SYYJ-SYYA,'b');
plot(0:1:Nt,SYYJRK-SYYA,'.g');
% plot(0:1:NtSYY3D-SYYA,'.k');
title('SIGMAyy')
% SIGMAxy'
subplot(3,1,3)
hold on
% plot(0:1:Nt,SXYJ-SXYA,'b');
% plot(0:1:Nt,SXYJRK-SXYA,'g');
plot(0:1:Nt,SXY3D-SXYA,'.k');
title('SIGMAxy')

