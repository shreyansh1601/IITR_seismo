% Solution of 2D Lagrangian temperature equation 
% on a regular grid with implicit finite-differences 
% for a moving medium with variable conductivity; 
% using of marker-in-cell approach for advection of temperature  

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
% Volumetric heat capacity, J/K/m^3
MRHOCP(1)=3200*1000; % background medium
MRHOCP(2)=3300*1100; % hot block
% thermal conductivity, W/m/K
MKT(1)=3; % background medium
MKT(2)=10; % hot block

% Subgrid diffusion coefficient
dsubgridt=1;

% Making vectors for nodal points positions (basic nodes)
gridx=0:xstp:xsize; % Horizontal
gridy=0:ystp:ysize; % Vertical

% max thermal diffusivity, m^2/s
kappa=max(MKT)/min(MRHOCP); 
% Timestep
dtexp=min(xstp,ystp)^2/3/kappa; % Limitation for explicit timestep
% Timestep limitation for advection
if(vx~=0)
    dtexp=min(dtexp,abs(xstp/vx)); % Limitation for horizontal advection timestep
end
if(vy~=0)
    dtexp=min(dtexp,abs(ystp/vy)); % Limitation for vertical advection timestep
end
timestep=0.5*dtexp; % Timestep, s

% Defining number of markers and steps between them in the horizontal and vertical direction
mxnum=200; %total number of markers in horizontal direction
mynum=300;  %total number of markers in vertical direction
mxstep=xsize/mxnum; %step between markers in horizontal direction   
mystep=ysize/mynum; %step between markers in vertical direction

% Creating markers arrays
MX=zeros(mynum*mxnum,1);   % X coordinate, m
MY=zeros(mynum*mxnum,1);   % Y coordinate, m
MI=zeros(mynum*mxnum,1);   % Type
MTK=zeros(mynum*mxnum,1);   % Temperature, K

% Defining intial position of markers
% Defining lithological structure of the model
% Marker counter
mm1=0;
for xm = 1:1:mxnum
    for ym = 1:1:mynum
        
        % Update marker counter:
        mm1=mm1+1;
        
        % Coordinates with small random displacement
        MX(mm1)=xm*mxstep-mxstep/2+(rand-0.5)*mxstep;
        MY(mm1)=ym*mystep-mystep/2+(rand-0.5)*mystep;
        
        if(MY(mm1)>ysize*0.3 && MY(mm1)<ysize*0.7 && MX(mm1)>xsize*0.3 && MX(mm1)<xsize*0.7) 
            % 2 = Hot block
            MI(mm1)=2;
            MTK(mm1)=twave;
        else
            % 1 = Medium
            MI(mm1)=1;
            MTK(mm1)=tback;
        end
    end
end

% Save Number of markers
marknum=mm1

% Density, heat capacity etc arrays
tk1 = zeros(ynum,xnum);     % Temperature
kt = zeros(ynum,xnum);     % thermal conductivity
rhocp = zeros(ynum,xnum);  % volumetric heat capacity

% Matrix of coefficients initialization for implicit solving
L=sparse(xnum*ynum,xnum*ynum);
% Vector of right part initialization for implicit solving
R=zeros(xnum*ynum,1);


% Time cycle
timesum=0; % Elapsed time
for ntimestep=1:1:tnum
    
    % Backup transport properties arrays
    rhocp0=rhocp; % volumetric heat capacity
    kt0=kt; % thermal conductivity 
    tk0=tk1; % temperature 
    
    % Clear transport properties arrays
    rhocp = zeros(ynum,xnum);    % volumetric heat capacity
    tk1 = zeros(ynum,xnum);     % Temperature
    kt = zeros(ynum,xnum);     % thermal conductivity
    % Clear wights for basic nodes
    wtnodes=zeros(ynum,xnum);

    % Interpolating parameters from markers to nodes
    for mm1 = 1:1:marknum

        % Check markers inside the grid
        if (MX(mm1)>=gridx(1) && MX(mm1)<=gridx(xnum) && MY(mm1)>=gridy(1) && MY(mm1)<=gridy(ynum)) 

            %  xn    rho(xn,yn)--------------------rho(xn+1,yn)
            %           ?           ^                  ?
            %           ?           ?                  ?
            %           ?          dy                  ?
            %           ?           ?                  ?
            %           ?           v                  ?
            %           ?<----dx--->o Mrho(xm,ym)       ?
            %           ?                              ?
            %           ?                              ?
            %  xn+1  rho(xn,yn+1)-------------------rho(xn+1,yn+1)
            %
            %
            % Define indexes for upper left node in the cell where the marker is
            xn=fix(MX(mm1)/xstp)+1;
            yn=fix(MY(mm1)/ystp)+1;
            if (xn<1)
                xn=1;
            end
            if (xn>xnum-1)
                xn=xnum-1;
            end
            if (yn<1)
                yn=1;
            end
            if (yn>ynum-1)
                yn=ynum-1;
            end

            % Define normalized distances from marker to the upper left node;
            dx=(MX(mm1)-gridx(xn))./xstp;
            dy=(MY(mm1)-gridy(yn))./ystp;
            % Compute weights for 4 surrounding nodes
            wtmij=(1-dx)*(1-dy);
            wtmi1j=(1-dx)*(dy);    
            wtmij1=(dx)*(1-dy);
            wtmi1j1=(dx)*(dy);

            % Define material volumetric heat capacity, etc from marker type
            MRHOCPCUR=MRHOCP(MI(mm1)); % volumetric heat capacity
            MKTCUR=MKT(MI(mm1)); % thermal conductivity

            % Add properties to 4 surrounding basic nodes
            % Upper-Left node
            rhocp(yn,xn)=rhocp(yn,xn)+wtmij*MRHOCPCUR;
            kt(yn,xn)=kt(yn,xn)+wtmij*MKTCUR;
            tk1(yn,xn)=tk1(yn,xn)+wtmij*MRHOCPCUR*MTK(mm1);
            wtnodes(yn,xn)=wtnodes(yn,xn)+wtmij;
            % Lower-Left node
            rhocp(yn+1,xn)=rhocp(yn+1,xn)+wtmi1j*MRHOCPCUR;
            kt(yn+1,xn)=kt(yn+1,xn)+wtmi1j*MKTCUR;
            tk1(yn+1,xn)=tk1(yn+1,xn)+wtmi1j*MRHOCPCUR*MTK(mm1);
            wtnodes(yn+1,xn)=wtnodes(yn+1,xn)+wtmi1j;
            % Upper-Right node
            rhocp(yn,xn+1)=rhocp(yn,xn+1)+wtmij1*MRHOCPCUR;
            kt(yn,xn+1)=kt(yn,xn+1)+wtmij1*MKTCUR;
            tk1(yn,xn+1)=tk1(yn,xn+1)+wtmij1*MRHOCPCUR*MTK(mm1);
            wtnodes(yn,xn+1)=wtnodes(yn,xn+1)+wtmij1;
            % Lower-Right node
            rhocp(yn+1,xn+1)=rhocp(yn+1,xn+1)+wtmi1j1*MRHOCPCUR;
            kt(yn+1,xn+1)=kt(yn+1,xn+1)+wtmi1j1*MKTCUR;
            tk1(yn+1,xn+1)=tk1(yn+1,xn+1)+wtmi1j1*MRHOCPCUR*MTK(mm1);
            wtnodes(yn+1,xn+1)=wtnodes(yn+1,xn+1)+wtmi1j1;

        end
    end
    
    % Computing  density, thermal conductivity etc for nodal points
    for i=1:1:ynum
        for j=1:1:xnum
            % Density
            if (wtnodes(i,j)>0)
                % Compute new value interpolated from markers
                tk1(i,j)=tk1(i,j)./rhocp(i,j);
                rhocp(i,j)=rhocp(i,j)./wtnodes(i,j);
                kt(i,j)=kt(i,j)./wtnodes(i,j);
            else
                % If no new value is interpolated from markers old value is used
                rhocp(i,j)=rhocp0(i,j);
                kt(i,j)=kt0(i,j);
                tk1(i,j)=tk0(i,j);
            end
        end
    end
    
    % Applying thermal boundary conditions for interpolated temperature
    % Upper boundary 
%     tk1(1,:)=tback; % Constant temperature
    tk1(1,:)=tk1(2,:); % Insulating boundary
    % Lower boundary 
%     tk1(ynum,:)=tback; % Constant temperature
    tk1(ynum,:)=tk1(ynum-1,:); % Insulating boundary
    % Left boundary
%     tk1(2:ynum-1,1)=tback; % Constant temperature
    tk1(2:ynum-1,1)=tk1(2:ynum-1,2); % Insulating boundary
    % Right boundary
%     tk1(2:ynum-1,xnum)=tback; % Constant temperature
    tk1(2:ynum-1,xnum)=tk1(2:ynum-1,xnum-1); % Insulating boundary
     

    
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
%                     % Constant temperature: T(i,j)=tback
%                     L(k,k)=1;
%                     R(k,1)=tback;
                    % Insulating boundary: T(i,j)=T(i+1,j)
                    L(k,k)=1;
                    L(k,k+1)=-1;
                    R(k,1)=0;                    
                end
                % Lower boundary
                if(i==ynum)
%                     % Constant temperature: T(i,j)=tback
%                     L(k,k)=1;
%                     R(k,1)=tback;
                    % Insulating boundary: T(i,j)=T(i-1,j)
                    L(k,k)=1;
                    L(k,k-1)=-1;
                    R(k,1)=0;                    
                end
                % Left boundary
                if(j==1 && i>1 && i<ynum)
%                     % Constant temperature: T(i,j)=tback
%                     L(k,k)=1;
%                     R(k,1)=tback;
                    % Insulating boundary: T(i,j)=T(i,j+1)
                    L(k,k)=1;
                    L(k,k+ynum)=-1;
                    R(k,1)=0;                    
                end
                % Right boundary
                if(j==xnum && i>1 && i<ynum)
%                     % Constant temperature: T(i,j)=tback
%                     L(k,k)=1;
%                     R(k,1)=tback;
                    % Insulating boundary: T(i,j)=T(i,j-1)
                    L(k,k)=1;
                    L(k,k-ynum)=-1;
                    R(k,1)=0;                    
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
                L(k,k)=L(k,k)+rhocp(i,j)/timestep; % ADD coefficient for T(i,j)
                % Right part
                R(k,1)=rhocp(i,j)*tk1(i,j)/timestep;
            end
        end
    end
    % Obtaining solution
    S=L\R;
    % Reloading solution to the new temperature array
    % Reloading solution to the new temperature array
    tk2=tk1;
    for i=1:1:ynum
        for j=1:1:xnum
            % Global index
            k=(j-1)*ynum+i;
            % Reload solution
            tk2(i,j)=S(k);
        end
    end
    % Computing temperature changes
    dtk1=tk2-tk1;
    
    % Computing subgrid diffusion for markers
    if (dsubgridt>0)
        % Clear subgrid temperature changes for nodes
        dtkn=zeros(ynum,xnum);
        % Clear wights for basic nodes
        wtnodes=zeros(ynum,xnum);
        % Marker cycle
        for mm1 = 1:1:marknum
            % Check markers inside the grid
            if (MX(mm1)>=0 && MX(mm1)<=xsize && MY(mm1)>=0 && MY(mm1)<=ysize) 

                %  yn    T(yn,xn)--------------------T(yn,xn+1)
                %           ?           ^                  ?
                %           ?           ?                  ?
                %           ?          dy                  ?
                %           ?           ?                  ?
                %           ?           v                  ?
                %           ?<----dx--->o MTK(mm1)         ?
                %           ?                              ?
                %           ?                              ?
                %  yn+1  T(yn+1,xn)-------------------V(yn+1,xn+1)
                %
                %
                % Interpolating temperature changes from basic nodes
                xn=fix(MX(mm1)/xstp)+1;
                yn=fix(MY(mm1)/ystp)+1;
                % Check horizontal index for upper left T-node 
                % It must be between 1 and xnum-1 (see picture for grid)
                if (xn<1)
                    xn=1;
                end
                if (xn>xnum-1)
                    xn=xnum-1;
                end
                if (yn<1)
                    yn=1;
                end
                if (yn>ynum-1)
                    yn=ynum-1;
                end
                % Define and check normalized distances from marker to the upper left T-node;
                dx=MX(mm1)./xstp-xn+1;
                dy=MY(mm1)./ystp-yn+1;
                % Compute weights for 4 surrounding nodes
                wtmij=(1-dx)*(1-dy);
                wtmi1j=(1-dx)*(dy);    
                wtmij1=(dx)*(1-dy);
                wtmi1j1=(dx)*(dy);
                % Interpolate old nodal temperature for the marker
                tkm=wtmij*tk1(yn,xn)+wtmi1j*tk1(yn+1,xn)+...
                    wtmij1*tk1(yn,xn+1)+wtmi1j1*tk1(yn+1,xn+1);
                % Calculate Nodal-Marker subgrid temperature difference
                dtkm=tkm-MTK(mm1);

                % Compute local thermal diffusion timescale for the marker
                MRHOCPCUR=MRHOCP(MI(mm1)); % volumetric heat capacity
                MKTCUR=MKT(MI(mm1)); % thermal conductivity
                tdm=MRHOCPCUR/MKTCUR/(2/xstp^2+2/ystp^2);

                % Computing subgrid diffusion
                sdif=-dsubgridt*timestep/tdm;
                if(sdif<-30) 
                    sdif=-30;
                end
                dtkm=dtkm*(1-exp(sdif));    

                % Correcting old temperature for the marker
                MTK(mm1)=MTK(mm1)+dtkm;

                % Interpolating subgrid temperature changes 
                % to 4 surrounding nodes
                dtkn(yn,xn)=dtkn(yn,xn)+MRHOCPCUR*wtmij*dtkm;
                wtnodes(yn,xn)=wtnodes(yn,xn)+MRHOCPCUR*wtmij;
                dtkn(yn+1,xn)=dtkn(yn+1,xn)+MRHOCPCUR*wtmi1j*dtkm;
                wtnodes(yn+1,xn)=wtnodes(yn+1,xn)+MRHOCPCUR*wtmi1j;
                dtkn(yn,xn+1)=dtkn(yn,xn+1)+MRHOCPCUR*wtmij1*dtkm;
                wtnodes(yn,xn+1)=wtnodes(yn,xn+1)+MRHOCPCUR*wtmij1;
                dtkn(yn+1,xn+1)=dtkn(yn+1,xn+1)+MRHOCPCUR*wtmi1j1*dtkm;
                wtnodes(yn+1,xn+1)=wtnodes(yn+1,xn+1)+MRHOCPCUR*wtmi1j1;
            end
        end

        % Computing subgrid diffusion for nodes
        for i=1:1:ynum
            for j=1:1:xnum
                % Density
                if (wtnodes(i,j)>0)
                    % Compute new value interpolated from markers
                    dtkn(i,j)=dtkn(i,j)./wtnodes(i,j);
                end
            end
        end

        % Subtracting subgrid diffusion part from nodal temperature changes
        dtk1=dtk1-dtkn;

    end

    % Updating temperature for markers
    for mm1 = 1:1:marknum
        % Check markers inside the grid
        if (MX(mm1)>=0 && MX(mm1)<=xsize && MY(mm1)>=0 && MY(mm1)<=ysize) 

            %  yn    T(yn,xn)--------------------T(yn,xn+1)
            %           ?           ^                  ?
            %           ?           ?                  ?
            %           ?          dy                  ?
            %           ?           ?                  ?
            %           ?           v                  ?
            %           ?<----dx--->o MTK(mm1)         ?
            %           ?                              ?
            %           ?                              ?
            %  yn+1  T(yn+1,xn)-------------------V(yn+1,xn+1)
            %
            %
            % Interpolating temperature changes from basic nodes
            xn=fix(MX(mm1)/xstp)+1;
            yn=fix(MY(mm1)/ystp)+1;
            % Check horizontal index for upper left T-node 
            % It must be between 1 and xnum-1 (see picture for grid)
            if (xn<1)
                xn=1;
            end
            if (xn>xnum-1)
                xn=xnum-1;
            end
            if (yn<1)
                yn=1;
            end
            if (yn>ynum-1)
                yn=ynum-1;
            end
            % Define and check normalized distances from marker to the upper left T-node;
            dx=MX(mm1)./xstp-xn+1;
            dy=MY(mm1)./ystp-yn+1;
            % Compute weights for 4 surrounding nodes
            wtmij=(1-dx)*(1-dy);
            wtmi1j=(1-dx)*(dy);    
            wtmij1=(dx)*(1-dy);
            wtmi1j1=(dx)*(dy);
            % Calculate Marker temperature change from four surrounding nodes
            dtkm=wtmij*dtk1(yn,xn)+wtmi1j*dtk1(yn+1,xn)+...
            	wtmij1*dtk1(yn,xn+1)+wtmi1j1*dtk1(yn+1,xn+1);
            %
            % Computing new temperature for the marker
            if(ntimestep>1)
                % Add temperature change
                MTK(mm1)=MTK(mm1)+dtkm;
            else
                % Interpolate new temperature for the 1st timestep
                MTK(mm1)=wtmij*tk2(yn,xn)+wtmi1j*tk2(yn+1,xn)+...
                        wtmij1*tk2(yn,xn+1)+wtmi1j1*tk2(yn+1,xn+1);
            end

        end
    end
   
    
    
    % Moving Markers by velocity field
    % Create arrays for velocity of markers
    vxm=zeros(4,1);
    vym=zeros(4,1);
    % Marker cycle
    for mm1 = 1:1:marknum

        % Displacing Marker according to its velocity
        MX(mm1)=MX(mm1)+timestep*vx;
        MY(mm1)=MY(mm1)+timestep*vy;
        % Recycle markers
        if(MX(mm1)<0)
            MX(mm1)=MX(mm1)+xsize;
        end
        if(MX(mm1)>xsize)
            MX(mm1)=MX(mm1)-xsize;
        end
        if(MY(mm1)<0)
            MY(mm1)=MY(mm1)+ysize;
        end
        if(MY(mm1)>ysize)
            MY(mm1)=MY(mm1)-ysize;
        end

    end
    
    
    % Open figure
    figure(1);
    % Plotting implicit solution
    subplot(1,2,1);
    pcolor(gridx/1000,gridy/1000,tk2);
    shading interp;
    colorbar;
    axis ij image;
    title(['Temperature, K: step=',num2str(ntimestep),' time,Myr=',num2str(timesum/(1e+6*365.25*24*3600))]);
    % Plotting density
    subplot(1,2,2);
    pcolor(gridx/1000,gridy/1000,rhocp);
    shading interp;
    colorbar;
    axis ij image;
    title(['RHO*CP, J/K/m^3: step=',num2str(ntimestep),' time,Myr=',num2str(timesum/(1e+6*365.25*24*3600))]);
    
    % Stop for .1 second
    pause(.1);
    
    % Add time counter
    timesum=timesum+timestep;
end