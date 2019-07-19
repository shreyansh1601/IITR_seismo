% Comparison of upwind, downwind and central differences 
% for 1D advection of a square density wave in a constant velocity field

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
vx=1;
dt=xstp/vx*60/800;
ntimes=800;

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

% Arrays for upwind differences
densnewup=densold; 
densoldup=densold;
% Arrays for downwind differences
densnewdn=densold; 
densolddn=densold;
% Arrays for central differences
densnewcn=densold; 
densoldcn=densold;

% Open Figure
figure(1);
% Advect density with various schemes 
for t=1:1:ntimes
    for i=2:1:xnum-1
        if(vx>0)
            % Upwind
            densnewup(i)=densoldup(i)-vx*dt*(densoldup(i)-densoldup(i-1))/xstp;
            % Downwind
            densnewdn(i)=densolddn(i)-vx*dt*(densolddn(i+1)-densolddn(i))/xstp;
        else
            % Upwind
            densnewup(i)=densoldup(i)-vx*dt*(densoldup(i+1)-densoldup(i))/xstp;
            % Downwind
            densnewdn(i)=densolddn(i)-vx*dt*(densolddn(i)-densolddn(i-1))/xstp;
        end
        %Central
        densnewcn(i)=densoldcn(i)-vx*dt*(densoldcn(i+1)-densoldcn(i-1))/2/xstp;
    end
    % Reset old values for the next step
    densoldup=densnewup;
    densolddn=densnewdn;
    densoldcn=densnewcn;
    
    % Plot results 
    % Upwind
    subplot(1,3,1)
    plot(xgrid,densoldup);
    title('upwind');
    %defining horizontal and vertical axis
    axis([0 xsize 2000 4000]);
    % Downwind 
    subplot(1,3,2)
    plot(xgrid,densolddn);
    title('downwind');
    %defining horizontal and vertical axis
    axis([0 xsize 2000 4000]);
    % Central 
    subplot(1,3,3)
    plot(xgrid,densoldcn);
    title('cetral');
    %defining horizontal and vertical axis
    axis([0 xsize 2000 4000]);
    pause(0.01);
end


