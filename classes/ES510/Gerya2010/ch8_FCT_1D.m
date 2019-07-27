% Using of FCT algorithm (Boris and Book, 1973) 
% for 1D advection of a square density wave in constant velocity field

% Clear all arrays
clear all; 

% Clear all figures
clf;


% Defining model size and resolution 
xsize=150.0;
xnum=151;
xstp=xsize/(xnum-1);

% Defining numerical grid 
xgrid=0:xstp:xsize;

% Defining advection velocity and timestep
vx=.1;
dt=xstp/abs(vx)*60/800
ntimes=800;
eps=vx*dt/xstp;
nu=1/8+(eps^2)/2;
mu=1/8;

% Defining initial density distribution 
for i=1:1:xnum
    % Background density
    densold(i)=3000;
    % Square wave
    if (xgrid(i)>=3 && xgrid(i)<=23)
        densold(i)=3300;
    end
    % Triangular wave
    if (xgrid(i)>=43 && xgrid(i)<=53)
        densold(i)=3000+(xgrid(i)-43)/10*300;
    end
    if (xgrid(i)>=53 && xgrid(i)<=63)
        densold(i)=3300-(xgrid(i)-53)/10*300;
    end
end
densnew=densold;
denscor=densold;

% Open Figure
figure(1);
% Advect density
for t=1:1:ntimes
    % Step 1: Transport+numerical diffusion stage
    for i=2:1:xnum-1
        densnew(i)=densold(i)-eps/2*(densold(i+1)-densold(i-1))+nu*(densold(i+1)-2*densold(i)+densold(i-1));
    end
    % Step 2: Antidiffusion stage
    % Antidiffusion flow for the first cell
    fadc(1)=0;
    for i=2:1:xnum-2
        % Corrected antidiffusion flow for current cell
        delt0=densnew(i)-densnew(i-1);
        delt1=densnew(i+1)-densnew(i);
        delt2=densnew(i+2)-densnew(i+1);
        s=sign(delt1);
        fadc(i)=s*max(0,min(min(s*delt2,s*delt0),mu*abs(delt1)));
        denscor(i)=densnew(i)-fadc(i)+fadc(i-1);
    end
    
densold=denscor;
% Plot results            
plot(xgrid,densold);
%defining horizontal and vertical axis
axis([0 xsize 2950 3350]);
pause(0.01);
end


