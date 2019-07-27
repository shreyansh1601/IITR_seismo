% Using of marker-in-cell algorithm with regular Eulerian grid
% for 1D advection of a square density wave in a constant velocity field

% Clear all arrays
clear all; 

% Clear all figures
clf;

% EULERIAN GRID
% Defining model size 
xsize=150.0;
% Eulerian grid resolution 
xnum=151;
% Grid Step
xstp=xsize/(xnum-1);
% Nodal point positions 
xgrid=0:xstp:xsize;


% LAGRANGIAN MARKERS
% Defining number of markers per Eulerian cell
mcel=10;
% Total number of markers
mnum=(xnum-1)*mcel;
% Distance between markers
mxstp=xstp/mcel
% Marker points positions 
xmark=mxstp/2:mxstp:xsize-mxstp/2;


% Defining advection velocity and timestep
vx=1;
dt=xstp/abs(vx)*60/800;
ntimes=800;

% Defining initial density distribution for markers
for i=1:1:mnum
    % Background density
    dmark(i)=3000;
    % Square wave
    if (xmark(i)>=3 && xmark(i)<=23)
        dmark(i)=3300;
    end
    % Triangular wave
    if (xmark(i)>=43 && xmark(i)<=53)
        dmark(i)=3000+(xmark(i)-43)/10*300;
    end
    if (xmark(i)>=53 && xmark(i)<=63)
        dmark(i)=3300-(xmark(i)-53)/10*300;
    end
end

% Open Figure
figure(1);
% Advect density with markers
for t=1:1:ntimes
    % Recomputing coordinates of markers
    for i=1:1:mnum
        xmark(i)=xmark(i)+vx*dt;
        % Recycle markers exiting model
        % Exit from the left -> recycle to the right
        if (xmark(i)<0) 
            xmark(i)=xsize+xmark(i);
        end
        % Exit from the right -> recycle to the left
        if (xmark(i)>xsize) 
            xmark(i)=xmark(i)-xsize;
        end
    end

% Recompute density on nodes
% Reset density and weights on nodes 
dnodes0=zeros(xnum);
wtnodes=zeros(xnum);
% Interpolating density from markers to nodes
for i=1:1:mnum
    % Check markers inside the grid
    if (xmark(i)>=0 && xmark(i)<=xsize)
        % Define index for nearest node from the left of current marker
        % Note that 0.5 is subtracted since int16(0.5)=1 int16(0.49)=0
        xn=double(int16(xmark(i)./xstp-0.5))+1;
        % Check index
        if (xn<1)
             xn=1;
        end
        if (xn>xnum-1)
            xn=xnum-1;
        end
        % Compute relative distances from the marker to the node
        dx=xmark(i)./xstp-xn+1;
        % Add Marker density to the nearest of two surrounding nodes
        % -+--------[xn]-----(i)----+--------[xn+1]--------+--------
        if(dx<0.5)
            dnodes0(xn)=dnodes0(xn)+(1.0-dx).*dmark(i);
            wtnodes(xn)=wtnodes(xn)+(1.0-dx);
        else
            dnodes0(xn+1)=dnodes0(xn+1)+dx.*dmark(i);
            wtnodes(xn+1)=wtnodes(xn+1)+dx;
        end
    end
end
% Recomputing density on nodes
for i=1:1:xnum
    % Old density remains on nodes 
    % where no new density is interpolated from markers
    if (wtnodes(i)>0)
        dnodes(i)=dnodes0(i)/wtnodes(i);
    end
end

% Plot results            
plot(xgrid,dnodes);
%plot(xmark,dmark);
%defining horizontal and vertical axis
axis([0 xsize 2950 3350]);
pause(0.01);
end


