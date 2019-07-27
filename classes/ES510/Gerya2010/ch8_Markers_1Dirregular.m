% Using of marker-in-cell algorithm with irregular Eulerian grid
% for 1D advection of a square density wave in a variable velocity field;
% using of bisection algorithm

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
% Randomise slightly positions of internal nodes
for i=2:1:xnum-1
    xgrid(i)=xgrid(i)+(rand-0.5)*xstp/4;
end

% LAGRANGIAN MARKERS
% Defining number of markers per Eulerian cell
mcel=10;
% Total number of markers
mnum=(xnum-1)*mcel;
% Distance between markers
mxstp=xstp/mcel;
% Marker points positions 
xmark=mxstp/2:mxstp:xsize-mxstp/2;
% Randomise slightly positions of markers
for i=1:1:mnum
    xmark(i)=xmark(i)+(rand-0.5)*mxstp/3;
end

% Defining advection velocity and timestep
vx=1;
dt=xstp/abs(vx)*60/800;
ntimes=800;

% Prescrie slightly randomised velocity on nodes
for i=1:1:xnum
    vxg(i)=vx*(1+(rand-0.5));
end

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
        % Define index for the nearest node to the left by bisection
        xnmin=1;
        xnmax=xnum;
        while ((xnmax-xnmin)>1)
            % !!! SUBTRACT 0.5 since int16(0.5)=1
            xn=double(int16((xnmax+xnmin)./2-0.5));
            if(xgrid(xn)>xmark(i))
                xnmax=xn;
            else
                xnmin=xn;
            end
        end
        xn=xnmin;
        % Check horizontal index
        if (xn<1)
            xn=1;
        end
        if (xn>xnum-1)
            xn=xnum-1;
        end
        % Define normalized distance from marker to the node
        dx=(xmark(i)-xgrid(xn))/(xgrid(xn+1)-xgrid(xn));
        % Define velocity for the marker from two surrounding nodes
        vxm=vxg(xn)*(1-dx)+vxg(xn+1)*dx;
        % Move marker
        xmark(i)=xmark(i)+vxm*dt;
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
        % Define index for the nearest node to the left by bisection
        xnmin=1;
        xnmax=xnum;
        while ((xnmax-xnmin)>1)
            % !!! SUBTRACT 0.5 since int16(0.5)=1
            xn=double(int16((xnmax+xnmin)./2-0.5));
            if(xgrid(xn)>xmark(i))
                xnmax=xn;
            else
                xnmin=xn;
            end
        end
        xn=xnmin;
        % Check horizontal index
        if (xn<1)
            xn=1;
        end
        if (xn>xnum-1)
            xn=xnum-1;
        end
        % Define normalized distance from marker to the node
        dx=(xmark(i)-xgrid(xn))/(xgrid(xn+1)-xgrid(xn));
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


