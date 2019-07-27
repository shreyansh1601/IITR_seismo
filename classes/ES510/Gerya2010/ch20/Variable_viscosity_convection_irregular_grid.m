% Solution of Stokes and continuity equations 
% with variable viscosity in 2D with direct solver
% by using external function Stokes_Continuity_solver_channel()
% solving temperature equation
% Setup corresponds to constant viscosity channel test
% with non steady temperature distribution
% 
% Staggered Grid 
% 
%     vx       vx       vx    
%
% vy  T---vy---T---vy---T   vy
%     |        |        |
%     vx   P   vx   P   vx    
%     |        |        |
% vy  T---vy---T---vy---T   vy
%     |        |        |
%     vx   P   vx   P   vx    
%     |        |        |
% vy  T---vy---T---vy---T   vy
%
%     vx       vx       vx    
% 
% Lines show basic grid
% Ghost nodes shown outside the basic grid
% are used for boundary conditions

% Clearing all variables and arrays
clear all;
% Clearing figures
clf;

% Temperature at the top, and bottom of the model, K
ttop=273;
tbottom=1273;


% Acceleration of Gravity, m/s^2
gx=0;
gy=10;

% Gas constant J/mol/K
RGAS=8.314;

% Model size, m
xsize=1000000;
ysize=1000000;

% Defining resolution
xnum=51;
ynum=51;

% Convecting medium
% Rock density (kg/m3): RHO*exp(-ALP*(T-T0))
MRHO(1,1)=4000;               % standard density, kg/m^3
MRHO(1,2)=2.5e-5;             % thermal expansion, 1/K
MRHO(1,3)=ttop;   % T0
% Rock flow law: 
% 0 = constant viscosity;
% 1 = isothermal power-law 2*EPSILONii=C1*SIGMAii^n
% 2 = temperature dependent ETA=A*exp(Ea/RT0(1-(T-T0)/T0))
% 3 = temperature & depth dependent ETA=A*exp(-b*(T-ttop)/(tbottom-ttop) + c*Y/yzize)
MLAW(1)=3;
META(1,1)=1e+23;     % A, Pa s
META(1,2)=0; % b
META(1,3)=0;         % c
% Rock heat capacity, J/K/kg
MCP(1)=1250;
% Rock thermal conductivity, W/m/K
MKT(1,1)=5;
% Maximal timestep, s
timemax=5e+5*(365.25*24*3600);
% Maximal marker displacement step, number of gridsteps
markmax=0.5;
% Moving Markers: 
% 0 = not moving at all
% 1 = simple advection
% 4 = 4-th order in space  Runge-Kutta
markmove=4;
% Velocity calculation
% 0 = by Solving momentum and continuity equations
% 1 = solid body rotation
movemod=0;
% Maximal temperature change, allowed for one timestep, K
tempmax=70;
% Amount of timesteps
stepmax=10000;

% Numerical Subgrid temperature diffusion coefficient
dsubgridt=1;
% Shear heating on(1)/off(0)
frictyn=0;

% Pressure boundary conditions
% prfirst(1) = boundary condition mode:
% 0 - pressure in one cell definition
% 1 - pressure at the top and in the bottom of the channel
prfirst(1)=0;
% prfirst(2) = boundary condition value
prfirst(2)=0;


% Velocity Boundary condition specified by bleft,bright,btop,bbot 
% are implemented from ghost nodes 
% directly into Stokes and continuity equations

% Upper, Lower boundaries: Free slip
for j=1:1:xnum+1
    % Upper boundary: 
    % vx(1,j)=btop(j,1)+vx(2,j)*btop(j,2)
    btop(j,1)=0;
    btop(j,2)=1;
    % vy(1,j)=btop(j,3)+vy(2,j)*btop(j,4)
    btop(j,3)=0;
    btop(j,4)=0;
    % Lower boundary:  
    % vx(ynum+1,j)=bbottom(j,1)+vx(ynum,j)*bbottom(j,2)
    bbottom(j,1)=0;
    bbottom(j,2)=1;
    % vy(ynum,j)=bbottom(j,3)+vy(ynum-1,j)*bbottom(j,4)
    bbottom(j,3)=0;
    bbottom(j,4)=0;
end

% Left, Right boundaries: Free slip
for i=1:1:ynum+1
    % Left boundary: 
    % vx(i,1)=bleft(i,1)+vx(i,2)*bleft(i,2)
    bleft(i,1)=0;
    bleft(i,2)=0;
    % vy(i,1)=bleft(i,3)+vy(i,2)*bleft(i,42)
    bleft(i,3)=0;
    bleft(i,4)=1;
    % Right boundary: 
    % vx(i,xnum)=bright(i,1)+vx(i,xnum-1)*bbright(i,2)
    bright(i,1)=0;
    bright(i,2)=0;
    % vy(i,xnum+1)=bright(i,3)+vx(i,xnum)*bbright(i,4)
    bright(i,3)=0;
    bright(i,4)=1;
end

% Defining average gridsteps
xstp=xsize./(xnum-1);
ystp=ysize./(ynum-1);

% % Defining gridline positions for regular basic grid
% % In horizontal direction
% gridx=zeros(xnum,1);
% for i=2:1:xnum
%     gridx(i)=gridx(i-1)+xstp;
% end
% % In vertical direction
% gridy=zeros(ynum,1);
% for i=2:1:ynum
%     gridy(i)=gridy(i-1)+ystp;
% end


% Defining gridline positions for 51x51 irregular basic grid
% In horizontal direction
gridy=zeros(ynum,1);
for i=2:1:11
    gridy(i)=gridy(i-1)+10000;
end
for i=12:1:16
    gridy(i)=gridy(i-1)+10000*1.2407234^(i-11);
end
for i=17:1:36
    gridy(i)=gridy(i-1)+30000;
end
for i=37:1:41
    gridy(i)=gridy(i-1)+10000*1.2407234^(42-i);
end
for i=42:1:51
    gridy(i)=gridy(i-1)+10000;
end
gridy(xnum)=ysize;
% In vertical direction
gridx=gridy*xsize/ysize;





% Thermal boundary conditions
% Upper, Lower boundaries: constant temperature
for j=1:1:xnum
    % Upper boundary
    % tk(1,j)=btopt(j,1)+tk(2,j)*btop(j,2)
    btopt(j,1)=ttop;
    btopt(j,2)=0;
    % Lower boundary
    % tk(ynum,j)=bbottomt(j,1)+tk(ynum-1,j)*bbottomt(j,2)
    bbottomt(j,1)=tbottom;
    bbottomt(j,2)=0;
end
% Left, Right boundaries: symmetry
for i=1:1:ynum
    % Left boundary
    % tk(i,1)=bleftt(i,1)+bleftt(i,2)*tk(i,2);
    bleftt(i,1)=0;
    bleftt(i,2)=1;
    % Right boundary
    % tk(i,xnum)=brightt(i,1)+brightt(i,2)*tk(i,xnum-1);
    brightt(i,1)=0;
    brightt(i,2)=1;
end


% Defining number of markers and steps between them in the horizontal and vertical direction
mxnum=200; %total number of markers in horizontal direction
mynum=200; %total number of markers in vertical direction
mxstep=xsize/mxnum; %step between markers in horizontal direction   
mystep=ysize/mynum; %step between markers in vertical direction

% Creating markers arrays
MX=zeros(mynum,mxnum);   % X coordinate, m
MY=zeros(mynum,mxnum);   % Y coordinate, m
MTK=zeros(mynum,mxnum);  % Temperature, K
MEII=zeros(mynum,mxnum); % Second invariant of strain rate, 1/s
MI=zeros(mynum,mxnum);   % Type
MXN=zeros(mynum,mxnum);  % Horizontal index
MYN=zeros(mynum,mxnum);  % Vertical index

% Defining intial position of markers
% Defining lithological structure of the model
for xm = 1:1:mxnum
    for ym = 1:1:mynum
        % Coordinates with small random displacement
        MX(ym,xm)=xm*mxstep-mxstep/2+(rand-0.5)*mxstep;
        MY(ym,xm)=ym*mystep-mystep/2+(rand-0.5)*mystep;
        % Material Type
        MI(ym,xm)=1;
%         % Initial temperature profile with lateral perturbation
%         xwt=MX(ym,xm)/xsize;
%         ywt=MY(ym,xm)/ysize;
%         %ymid=0.28;
%         ymid=0.5;
%         if(ywt>ymid)
%             ywt=(1-ymid)+ymid*(abs(ywt-ymid)/(1-ymid))^10;
%         else
%             ywt=(1-ymid)-(1-ymid)*(abs(ywt-ymid)/ymid)^10;
%         end
%         MTK(ym,xm)=(ttop+(tbottom-ttop)*ywt^(2+xwt)+tbottom+(ttop-tbottom)*(1-ywt)^(3-xwt))/2;
    end
end


% Rock type, density, viscosity, temperature, thermal conductivity and RHO*Cp arrays
typ1 = zeros(ynum,xnum);
etas1 = zeros(ynum,xnum);
etan1 = zeros(ynum-1,xnum-1);
rho1 = zeros(ynum,xnum);
tk1 = zeros(ynum,xnum);
tk2=tk1;
rhocp1 = zeros(ynum,xnum);
kt1 = zeros(ynum,xnum);

% Read temperature structure from data file
% fdata=fopen('data_1a_regular.txt','rt'); 
% fdata=fopen('data_1c_regular.txt','rt'); 
% fdata=fopen('data_2a_regular.txt','rt'); 
fdata=fopen('data_1a_irregular.txt','rt'); 
% fdata=fopen('data_1c_irregular.txt','rt'); 
% fdata=fopen('data_2a_irregular.txt','rt'); 
for i=1:1:51
    for j=1:1:51
        tkcur=fscanf(fdata,'%f',1);
        tk1(i,j)=tkcur;
    end
end
tk2=tk1;        
fclose(fdata);


% Initial time, s
timesum=0;

% Main Time cycle
for ntimestep=1:1:stepmax


    % Backup transport properties arrays
    typ0 = typ1;
    etas0 = etas1;
    etan0 = etan1;
    rho0 = rho1;
    tk0=tk2;
    rhocp0=rhocp1;
    kt0=kt1;
    % Clear transport properties arrays
    typ1 = zeros(ynum,xnum);
    etas1 = zeros(ynum,xnum);
    etan1 = zeros(ynum-1,xnum-1);
    rho1 = zeros(ynum,xnum);
    tk1 = zeros(ynum,xnum);
    rhocp1 = zeros(ynum,xnum);
    kt1 = zeros(ynum,xnum);
    % Clear wights for basic nodes
    wtnodes=zeros(ynum,xnum);
    % Clear wights for etas
    wtetas=zeros(ynum,xnum);
    % Clear wights for etan
    wtetan=zeros(ynum-1,xnum-1);

    % Computing grid steps for basic nodes
    xstp1=zeros(xnum-1,1);
    ystp1=zeros(ynum-1,1);
    for i=1:1:xnum-1
        xstp1(i)=gridx(i+1)-gridx(i);
    end
    for i=1:1:ynum-1
        ystp1(i)=gridy(i+1)-gridy(i);
    end  
    
    % Computing grids and grid steps for Vx, Vy nodes
    % Horizontal (for Vy)
    gridcx=zeros(xnum+1,1);
    xstpc1=zeros(xnum,1);
    % Vertical (for Vx)
    gridcy=zeros(ynum+1,1);
    ystpc1=zeros(ynum,1);
    % First and last nodes and steps (for external nodes)
    % Horizontal (for Vy)
    gridcx(1)=-xstp1(1)/2;
    xstpc1(1)=xstp1(1);
    gridcx(xnum+1)=xsize+xstp1(xnum-1)/2;
    xstpc1(xnum)=xstp1(xnum-1);
    % Vertical (for Vx)
    gridcy(1)=-ystp1(1)/2;
    ystpc1(1)=ystp1(1);
    gridcy(ynum+1)=ysize+ystp1(ynum-1)/2;
    ystpc1(ynum)=ystp1(ynum-1);
    % Internal nodes
    for i=2:1:xnum
        gridcx(i)=(gridx(i)+gridx(i-1))/2;
    end
    for i=2:1:ynum
        gridcy(i)=(gridy(i)+gridy(i-1))/2;
    end   
    % Internal grid steps
    for i=2:1:xnum-1
        xstpc1(i)=(gridx(i+1)-gridx(i-1))/2;
    end
    for i=2:1:ynum-1
        ystpc1(i)=(gridy(i+1)-gridy(i-1))/2;
    end   
    
    % Interpolating parameters from markers to nodes
    for xm = 1:1:mxnum
        for ym = 1:1:mynum 
            
            % Check markers inside the grid
            if (MX(ym,xm)>=0 && MX(ym,xm)<=xsize && MY(ym,xm)>=0 && MY(ym,xm)<=ysize) 

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
                % Define indexes for upper left node in the cell where the marker is
                % using bisection
                % Find horizontal index
                xnmin=1;
                xnmax=xnum;
                while ((xnmax-xnmin)>1)
                    % !!! SUBTRACT 0.5 since int16(0.5)=1
                    xn=double(int16((xnmax+xnmin)./2-0.5));
                    if(gridx(xn)>MX(ym,xm))
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
                % Save horizontal index
                MXN(ym,xm)=xn;
                
                % Find vertical index
                ynmin=1;
                ynmax=ynum;
                while ((ynmax-ynmin)>1)
                    % !!! SUBTRACT 0.5 since int16(0.5)=1
                    yn=double(int16((ynmax+ynmin)./2-0.5));
                    if(gridy(yn)>MY(ym,xm))
                        ynmax=yn;
                    else
                        ynmin=yn;
                    end
                end
                yn=ynmin;
                % Check vertical index
                if (yn<1)
                    yn=1;
                end
                if (yn>ynum-1)
                    yn=ynum-1;
                end
                % Save Vertical index
                MYN(ym,xm)=yn;

                % Define normalized distances from marker to the upper left node;
                dx=(MX(ym,xm)-gridx(xn))./xstp1(xn);
                dy=(MY(ym,xm)-gridy(yn))./ystp1(yn);
                
                % Compute marker weight koefficient from cell dimensions
                % Number of markers in a cell is in invert proportion to the cell volume
                mwt=1;%/xstp1(xn)/ystp1(yn);
                
                % Interpolate nodal temperature for the marker at time=0
                if(timesum==0)
                    tkm=0;
                    tkm=tkm+(1.0-dx).*(1.0-dy).*tk2(yn,xn);
                    tkm=tkm+(1.0-dx).*dy.*tk2(yn+1,xn);
                    tkm=tkm+dx.*(1.0-dy).*tk2(yn,xn+1);
                    tkm=tkm+dx.*dy.*tk2(yn+1,xn+1);
                    MTK(ym,xm)=tkm;
                end

                
                % Compute density from marker temperature
                MRHOCUR=MRHO(MI(ym,xm),1)*(1-MRHO(MI(ym,xm),2)*(MTK(ym,xm)-ttop));
 
                % Compute rho*Cp for marker 
                MRHOCPCUR=MRHOCUR*MCP(MI(ym,xm));

                % Compute thermal conductivity from marker temperature
                % Rock thermal conductivity, W/m/K
                MKTCUR=MKT(MI(ym,xm),1);
                
                % Add properties to 4 surrounding nodes
                rho1(yn,xn)=rho1(yn,xn)+(1.0-dx).*(1.0-dy).*MRHOCUR*mwt;
                tk1(yn,xn)=tk1(yn,xn)+(1.0-dx).*(1.0-dy).*MTK(ym,xm)*mwt;
                kt1(yn,xn)=kt1(yn,xn)+(1.0-dx).*(1.0-dy).*MKTCUR*mwt;
                rhocp1(yn,xn)=rhocp1(yn,xn)+(1.0-dx).*(1.0-dy).*MRHOCPCUR*mwt;
                wtnodes(yn,xn)=wtnodes(yn,xn)+(1.0-dx).*(1.0-dy)*mwt;
                
                rho1(yn+1,xn)=rho1(yn+1,xn)+(1.0-dx).*dy.*MRHOCUR*mwt;
                tk1(yn+1,xn)=tk1(yn+1,xn)+(1.0-dx).*dy.*MTK(ym,xm)*mwt;
                kt1(yn+1,xn)=kt1(yn+1,xn)+(1.0-dx).*dy.*MKTCUR*mwt;
                rhocp1(yn+1,xn)=rhocp1(yn+1,xn)+(1.0-dx).*dy.*MRHOCPCUR*mwt;
                wtnodes(yn+1,xn)=wtnodes(yn+1,xn)+(1.0-dx).*dy*mwt;
                
                rho1(yn,xn+1)=rho1(yn,xn+1)+dx.*(1.0-dy).*MRHOCUR*mwt;
                tk1(yn,xn+1)=tk1(yn,xn+1)+dx.*(1.0-dy).*MTK(ym,xm)*mwt;
                kt1(yn,xn+1)=kt1(yn,xn+1)+dx.*(1.0-dy).*MKTCUR*mwt;
                rhocp1(yn,xn+1)=rhocp1(yn,xn+1)+dx.*(1.0-dy).*MRHOCPCUR*mwt;
                wtnodes(yn,xn+1)=wtnodes(yn,xn+1)+dx.*(1.0-dy)*mwt;
                
                rho1(yn+1,xn+1)=rho1(yn+1,xn+1)+dx.*dy.*MRHOCUR*mwt;
                tk1(yn+1,xn+1)=tk1(yn+1,xn+1)+dx.*dy.*MTK(ym,xm)*mwt;
                kt1(yn+1,xn+1)=kt1(yn+1,xn+1)+dx.*dy.*MKTCUR*mwt;
                rhocp1(yn+1,xn+1)=rhocp1(yn+1,xn+1)+dx.*dy.*MRHOCPCUR*mwt;
                wtnodes(yn+1,xn+1)=wtnodes(yn+1,xn+1)+dx.*dy*mwt;

                % Computing Marker Viscosity
                % 0 = Constant Viscosity
                METACUR=META(MI(ym,xm),1);
                % 1 = Isothermal power-law 2*EPSILONii=C1*SIGMAii^n
                % ETA=(2/C1)^(1/n)/2/EPSii^(1-1/n)
                if (MLAW(MI(ym,xm))==1)
                    % ETA=(2/C1)^(1/n)/2/EPSii^(1-1/n)
                    C1=META(MI(ym,xm),1);
                    n=META(MI(ym,xm),2);
                    meiicur=MEII(ym,xm);
                    % Check second invariant of strain rate
                    if (meiicur<1e-20) 
                        meiicur=1e-20;
                    end
                    METACUR=C1^(-1/n)*(2*meiicur)^(1/n-1);
                end
                % 2 = temperature dependent ETA=A*exp(Ea/RT0(1-(T-T0)/T0))
                if (MLAW(MI(ym,xm))==2)
                    A=META(MI(ym,xm),1);
                    Ea=META(MI(ym,xm),2);
                    T0=META(MI(ym,xm),3);
                    T=MTK(ym,xm);
                    METACUR=A*exp(Ea/RGAS/T0*(1-(T-T0)/T0));
                end
                % 3 = temperature & depth dependent ETA=A*exp(-b*(T-ttop)/(tbottom-ttop) + c*Y/yzize)
                if (MLAW(MI(ym,xm))==3)
                    A=META(MI(ym,xm),1);
                    b=META(MI(ym,xm),2);
                    c=META(MI(ym,xm),3);
                    T=MTK(ym,xm);
                    Y=MY(ym,xm);
                    METACUR=A*exp(-b*(T-ttop)/(tbottom-ttop)+c*Y/ysize);
                end
                

                % Add shear viscosity etas() and rock type typ() to 4 surrounding nodes
                % only using markers located at <=0.5 gridstep distances from nodes
                if(dx<=0.5 && dy<=0.5)
                    etas1(yn,xn)=etas1(yn,xn)+(1.0-dx).*(1.0-dy).*METACUR*mwt;
                    typ1(yn,xn)=typ1(yn,xn)+(1.0-dx).*(1.0-dy).*MI(ym,xm)*mwt;
                    wtetas(yn,xn)=wtetas(yn,xn)+(1.0-dx).*(1.0-dy)*mwt;
                end
                if(dx<=0.5 && dy>=0.5)
                    etas1(yn+1,xn)=etas1(yn+1,xn)+(1.0-dx).*dy.*METACUR*mwt;
                    typ1(yn+1,xn)=typ1(yn+1,xn)+(1.0-dx).*dy.*MI(ym,xm)*mwt;
                    wtetas(yn+1,xn)=wtetas(yn+1,xn)+(1.0-dx).*dy*mwt;
                end
                if(dx>=0.5 && dy<=0.5)
                    etas1(yn,xn+1)=etas1(yn,xn+1)+dx.*(1.0-dy).*METACUR*mwt;
                    typ1(yn,xn+1)=typ1(yn,xn+1)+dx.*(1.0-dy).*MI(ym,xm)*mwt;
                    wtetas(yn,xn+1)=wtetas(yn,xn+1)+dx.*(1.0-dy)*mwt;
                end
                if(dx>=0.5 && dy>=0.5)
                    etas1(yn+1,xn+1)=etas1(yn+1,xn+1)+dx.*dy.*METACUR*mwt;
                    typ1(yn+1,xn+1)=typ1(yn+1,xn+1)+dx.*dy.*MI(ym,xm)*mwt;
                    wtetas(yn+1,xn+1)=wtetas(yn+1,xn+1)+dx.*dy*mwt;
                end

                % Add normal viscosity etan() to the center of current cell
                etan1(yn,xn)=etan1(yn,xn)+(1.0-abs(0.5-dx)).*(1.0-abs(0.5-dy)).*METACUR*mwt;
                wtetan(yn,xn)=wtetan(yn,xn)+(1.0-abs(0.5-dx)).*(1.0-abs(0.5-dy))*mwt;
            end

        end
    end

    % Computing  Viscosity, density, rock type for nodal points
    for i=1:1:ynum;
        for j=1:1:xnum;
            % Density
            if (wtnodes(i,j)~=0)
                % Compute new value interpolated from markers
                rho1(i,j)=rho1(i,j)./wtnodes(i,j);
                tk1(i,j)=tk1(i,j)./wtnodes(i,j);
                kt1(i,j)=kt1(i,j)./wtnodes(i,j);
                rhocp1(i,j)=rhocp1(i,j)./wtnodes(i,j);
            else
                % If no new value is interpolated from markers old value is used
                rho1(i,j)=rho0(i,j);
                tk1(i,j)=tk0(i,j);
                kt1(i,j)=kt0(i,j);
                rhocp1(i,j)=rhocp0(i,j);
            end
            % Shear viscosity and type
            if (wtetas(i,j)~=0)
                % Compute new value interpolated from markers
                etas1(i,j)=etas1(i,j)./wtetas(i,j);
                typ1(i,j)=typ1(i,j)./wtetas(i,j);
            else
                % If no new value is interpolated from markers old value is used
                etas1(i,j)=etas0(i,j);
                typ1(i,j)=typ0(i,j);
            end
            % Normal viscosity
            if (i<ynum && j<xnum)
                if (wtetan(i,j)~=0)
                    % Compute new value interpolated from markers
                    etan1(i,j)=etan1(i,j)./wtetan(i,j);
                else
                    % If no new value is interpolated from markers old value is used
                    etan1(i,j)=etan0(i,j);
                end
            end
        end
    end
    
    % Applying thermal boundary conditions for interpolated temperature
    % Upper, Lower boundaries
    for j=2:1:xnum-1
        % Upper boundary
        tk1(1,j)=btopt(j,1)+btopt(j,2)*tk1(2,j);
        % Lower boundary
        tk1(ynum,j)=bbottomt(j,1)+bbottomt(j,2)*tk1(ynum-1,j);
    end
    % Left, Right boundaries: constant temperature
    for i=1:1:ynum
        % Left boundary
        tk1(i,1)=bleftt(i,1)+bleftt(i,2)*tk1(i,2);
        % Right boundary
        tk1(i,xnum)=brightt(i,1)+brightt(i,2)*tk1(i,xnum-1);
    end


    % Computing right part of mechanical equations
    % x-Stokes
    RX1=zeros(ynum+1,xnum);
    % y-Stokes
    RY1=zeros(ynum,xnum+1);
    % continuity
    RC1=zeros(ynum-1,xnum-1);
    % Grid points cycle
    for i=1:1:ynum;
        for j=1:1:xnum;
            % Right part of x-Stokes Equation
            if(j>1 && i>1 && j<xnum)
                RX1(i,j)=-gx*(rho1(i,j)+rho1(i-1,j))/2;
            end
            % Right part of y-Stokes Equation
            if(j>1 && i>1 && i<ynum)
                RY1(i,j)=-gy*(rho1(i,j)+rho1(i,j-1))/2;
            end
        end
    end


    % Computing velocity field
    if (movemod==0)
        % Solving of Stokes and Continuity equations on nodes
        % and computing residuals
        % by calling function Stokes_Continuity_solver_ghost()
        [vx1,resx1,vy1,resy1,pr1,resc1]=Stokes_Continuity_solver_grid(prfirst,etas1,etan1,xnum,ynum,gridx,gridy,RX1,RY1,RC1,bleft,bright,btop,bbottom);
    end
    % Solid body rotation
    if (movemod==1)
        for i=1:1:ynum+1;
            for j=1:1:xnum+1;
                % Vx
                if(j<xnum+1)
                    % Relative distance of vx node from the model center
                    dx=((j-1)*xstp-xsize/2)/(xsize/2);
                    dy=((i-1.5)*ystp-ysize/2)/(xsize/2);
                    dr=(dx^2+dy^2)^0.5;
                    % Set vx
                    vx1(i,j)=-vyright*dy;
                end
                % Vy
                if(i<ynum+1)
                    % Relative distance of vy node from the model center
                    dx=((j-1.5)*xstp-xsize/2)/(xsize/2);
                    dy=((i-1)*ystp-ysize/2)/(xsize/2);
                    dr=(dx^2+dy^2)^0.5;
                    % Set vy
                    vy1(i,j)=vyright*dx;
                end
            end
        end
    end

    % Computing EPS'xx=-EPS'yy, EPSxy=EPSyx deviatoric strain rate tensor components from vx, vy
    exy = zeros(ynum,xnum);
    exx = zeros(ynum-1,xnum-1);
    % Grid points cycle
    for i=1:1:ynum;
        for j=1:1:xnum;
            % EPS'xx=-EPS'yy=1/2(dvx/dx-dvy/dy)
            if(i<ynum && j<xnum)
                exx(i,j)=0.5*((vx1(i+1,j+1)-vx1(i+1,j))/xstp1(j)-(vy1(i+1,j+1)-vy1(i,j+1))/ystp1(i));
            end
            % EPSxy=EPSyx=1/2(dvx/dy+dvy/dx)
            exy(i,j)=0.5*((vx1(i+1,j)-vx1(i,j))/ystpc1(i)+(vy1(i,j+1)-vy1(i,j))/xstpc1(j));
        end
    end

    % Defining displacement timestep
    timestep=timemax % initial displacement step
    % Check maximal velocity
    vxmax=max(abs(max(max(vx1))),abs(min(min(vx1))));
    vymax=max(abs(max(max(vy1))),abs(min(min(vy1))));
    % Check marker displacement step
    if (vxmax>0)
        if (timestep>markmax*xstp/vxmax);
            timestep=markmax*xstp/vxmax;
        end
    end
    if (vymax>0)
        if (timestep>markmax*ystp/vymax);
            timestep=markmax*ystp/vymax;
        end
    end
    timestep=timestep % final displacement step
    
    % Solving Temperature equation
    if (timestep>0 && tempmax>0)
        
        % Computing right part of temperature equation
        RT1=zeros(ynum,xnum);
        % Computing shear heating for Temperature nodes
        % Hs=4*etan*EPS'xx^2+4*etas*EPSxy^2
        % Grid points cycle
        for i=2:1:ynum-1;
            for j=2:1:xnum-1;
                % Shear heating on(1)/off(0)
                if(frictyn==1)
                    % Adding 4*etas*EPSxy^2
                    RT1(i,j)=RT1(i,j)+4*etas1(i,j)*exy(i,j)^2;
                    % Computing and adding 4*etan*EPS'xx^2
                    RT1(i,j)=RT1(i,j)+etan1(i-1,j-1)*exx(i-1,j-1)^2+etan1(i,j-1)*exx(i,j-1)^2+etan1(i-1,j)*exx(i-1,j)^2+etan1(i,j)*exx(i,j)^2;
                end
            end
        end

        % Solving Temperature equation with displacement timestep
        [tk2,rest]=Temperature_solver_grid(timestep,xnum,ynum,gridx,gridy,kt1,rhocp1,tk1,RT1,bleftt,brightt,btopt,bbottomt);
        % Computing temperature changes
        dtk1=tk2-tk1;
        % Checking temperature changes
        dtkmax=max(max(abs(dtk1)))
        % Repeating temperature solution if temperature changes are too big
        if(dtkmax>tempmax)
            % Computing reduced timestep
            timestep=timestep*tempmax/dtkmax;
            % Solving Temperature equation with reduced timestep
            [tk2,rest]=Temperature_solver_grid(timestep,xnum,ynum,gridx,gridy,kt1,rhocp1,tk1,RT1,bleftt,brightt,btopt,bbottomt);
            % Computing temperature changes
            dtk1=tk2-tk1;
        end

        % Computing subgrid diffusion for markers
        if (dsubgridt>0)
            % Clear subgrid temperature changes for nodes
            dtkn=zeros(ynum,xnum);
            % Clear wights for basic nodes
            wtnodes=zeros(ynum,xnum);
            % Marker cycle
            for xm = 1:1:mxnum
                for ym = 1:1:mynum  
                    
                % Check markers inside the grid
                    if (MX(ym,xm)>=0 && MX(ym,xm)<=xsize && MY(ym,xm)>=0 && MY(ym,xm)<=ysize) 
                        
                        %  yn    T(yn,xn)--------------------T(yn,xn+1)
                        %           ?           ^                  ?
                        %           ?           ?                  ?
                        %           ?          dy                  ?
                        %           ?           ?                  ?
                        %           ?           v                  ?
                        %           ?<----dx--->o Mrho(ym,xm)       ?
                        %           ?                              ?
                        %           ?                              ?
                        %  yn+1  T(yn+1,xn)-------------------V(yn+1,xn+1)
                        %
                        %
                        % Interpolating temperature changes from basic nodes
                        %
                        % Define indexes for upper left node in the cell where the marker is
                        xn=MXN(ym,xm);
                        yn=MYN(ym,xm);

                        % Define normalized distances from marker to the upper left node;
                        dx=(MX(ym,xm)-gridx(xn))./xstp1(xn);
                        dy=(MY(ym,xm)-gridy(yn))./ystp1(yn);

                        % Compute marker weight koefficient from cell dimensions
                        % Number of markers in a cell is in invert proportion to the cell volume
                        mwt=1;%/xstp1(xn)/ystp1(yn);


                        % Interpolate old nodal temperature for the marker
                        tkm=0;
                        tkm=tkm+(1.0-dx).*(1.0-dy).*tk1(yn,xn);
                        tkm=tkm+(1.0-dx).*dy.*tk1(yn+1,xn);
                        tkm=tkm+dx.*(1.0-dy).*tk1(yn,xn+1);
                        tkm=tkm+dx.*dy.*tk1(yn+1,xn+1);
                        % Calculate Nodal-Marker subgrid temperature difference
                        dtkm=tkm-MTK(ym,xm);
                        % Compute nodal k and RHO*Cp for the marker 
                        % k
                        ktm=0;
                        ktm=ktm+(1.0-dx).*(1.0-dy).*kt1(yn,xn);
                        ktm=ktm+(1.0-dx).*dy.*kt1(yn+1,xn);
                        ktm=ktm+dx.*(1.0-dy).*kt1(yn,xn+1);
                        ktm=ktm+dx.*dy.*kt1(yn+1,xn+1);
                        % RHO*Cp
                        rhocpm=0;
                        rhocpm=rhocpm+(1.0-dx).*(1.0-dy).*rhocp1(yn,xn);
                        rhocpm=rhocpm+(1.0-dx).*dy.*rhocp1(yn+1,xn);
                        rhocpm=rhocpm+dx.*(1.0-dy).*rhocp1(yn,xn+1);
                        rhocpm=rhocpm+dx.*dy.*rhocp1(yn+1,xn+1);

                        % Compute local thermal diffusion timescale for the marker
                        tdm=rhocpm/ktm/(2/xstp^2+2/ystp^2);

                        % Computing subgrid diffusion
                        sdif=-dsubgridt*timestep/tdm;
                        if(sdif<-30) 
                            sdif=-30;
                        end
                        dtkm=dtkm*(1-exp(sdif));    

                        % Correcting old temperature for the marker
                        MTK(ym,xm)=MTK(ym,xm)+dtkm;

                        % Interpolating subgrid temperature changes to 4 nodes
                        dtkn(yn,xn)=dtkn(yn,xn)+(1.0-dx).*(1.0-dy).*dtkm*mwt;
                        wtnodes(yn,xn)=wtnodes(yn,xn)+(1.0-dx).*(1.0-dy)*mwt;
                        
                        dtkn(yn+1,xn)=dtkn(yn+1,xn)+(1.0-dx).*dy.*dtkm*mwt;
                        wtnodes(yn+1,xn)=wtnodes(yn+1,xn)+(1.0-dx).*dy*mwt;
                        
                        dtkn(yn,xn+1)=dtkn(yn,xn+1)+dx.*(1.0-dy).*dtkm*mwt;
                        wtnodes(yn,xn+1)=wtnodes(yn,xn+1)+dx.*(1.0-dy)*mwt;
                        
                        dtkn(yn+1,xn+1)=dtkn(yn+1,xn+1)+dx.*dy.*dtkm*mwt;
                        wtnodes(yn+1,xn+1)=wtnodes(yn+1,xn+1)+dx.*dy*mwt;
                        
                    end
                end
            end

            % Computing subgrid diffusion for nodes
            for i=1:1:ynum;
                for j=1:1:xnum;
                    % Density
                    if (wtnodes(i,j)~=0)
                        % Compute new value interpolated from markers
                        dtkn(i,j)=dtkn(i,j)./wtnodes(i,j);
                    end
                end
            end
            
            % Subtracting subgrid diffusion part from nodal temperature changes
            dtk1=dtk1-dtkn;
            
        end
        
        % Updating temperature for markers
        for xm = 1:1:mxnum
            for ym = 1:1:mynum  

                % Check markers inside the grid
                if (MX(ym,xm)>=0 && MX(ym,xm)<=xsize && MY(ym,xm)>=0 && MY(ym,xm)<=ysize) 

                    %  yn    T(yn,xn)--------------------T(yn,xn+1)
                    %           ?           ^                  ?
                    %           ?           ?                  ?
                    %           ?          dy                  ?
                    %           ?           ?                  ?
                    %           ?           v                  ?
                    %           ?<----dx--->o Mrho(ym,xm)       ?
                    %           ?                              ?
                    %           ?                              ?
                    %  yn+1  T(yn+1,xn)-------------------V(yn+1,xn+1)
                    %
                    %
                    % Interpolating temperature changes from basic nodes
                    %
                    % Define indexes for upper left node in the cell where the marker is
                    xn=MXN(ym,xm);
                    yn=MYN(ym,xm);

                    % Define normalized distances from marker to the upper left node;
                    dx=(MX(ym,xm)-gridx(xn))./xstp1(xn);
                    dy=(MY(ym,xm)-gridy(yn))./ystp1(yn);

                    % Calculate Marker temperature change from four surrounding nodes
                    dtkm=0;
                    dtkm=dtkm+(1.0-dx).*(1.0-dy).*dtk1(yn,xn);
                    dtkm=dtkm+(1.0-dx).*dy.*dtk1(yn+1,xn);
                    dtkm=dtkm+dx.*(1.0-dy).*dtk1(yn,xn+1);
                    dtkm=dtkm+dx.*dy.*dtk1(yn+1,xn+1);
                    %
                    %Computing new temperature for the marker
                    MTK(ym,xm)=MTK(ym,xm)+dtkm;
                    
                end
            end
        end

    end
    
    

    % Computing strain rate invariant for markers
    for xm = 1:1:mxnum
        for ym = 1:1:mynum  

            % Check markers inside the grid
            if (MX(ym,xm)>=0 && MX(ym,xm)<=xsize && MY(ym,xm)>=0 && MY(ym,xm)<=ysize) 
            
                %  xn    V(xn,yn)--------------------V(xn+1,yn)
                %           ?           ^                  ?
                %           ?           ?                  ?
                %           ?          dy                  ?
                %           ?           ?                  ?
                %           ?           v                  ?
                %           ?<----dx--->o Mrho(xm,ym)       ?
                %           ?                              ?
                %           ?                              ?
                %  xn+1  V(xn,yn+1)-------------------V(xn+1,yn+1)

                % Define indexes for upper left BASIC node in the cell where the marker is
                xnmin=MXN(ym,xm);
                ynmin=MYN(ym,xm);
                
                % Calculating Second invariant of strain rate for marker
                %
                % Interpolating squares of EPS'xx=-EPS'yy from cell centers
                % EPS'xx-nodes are displaced rightward and downward for 1/2 of gridsteps
                % Horizontal EPS'xx index
                xn=xnmin;
                if(MX(ym,xm)<gridcx(xn+1))
                    xn=xn-1;
                end
                if (xn<1)
                    xn=1;
                end
                if (xn>xnum-2)
                    xn=xnum-2;
                end
                % Vertical EPS'xx index
                yn=ynmin;
                if(MY(ym,xm)<gridcy(yn+1))
                    yn=yn-1;
                end
                if (yn<1)
                    yn=1;
                end
                if (yn>ynum-2)
                    yn=ynum-2;
                end
               
                % Define and check normalized distances from marker to the upper left EPS'xx-node;
                dx=(MX(ym,xm)-gridcx(xn+1))./xstpc1(xn+1);
                dy=(MY(ym,xm)-gridcy(yn+1))./ystpc1(yn+1);
                
                % Calculate Marker EPS'xx^2 from four surrounding nodes
                exxm=0;
                exxm=exxm+(1.0-dx).*(1.0-dy).*exx(yn,xn)*exx(yn,xn);
                exxm=exxm+(1.0-dx).*dy.*exx(yn+1,xn)*exx(yn+1,xn);
                exxm=exxm+dx.*(1.0-dy).*exx(yn,xn+1)*exx(yn,xn+1);
                exxm=exxm+dx.*dy.*exx(yn+1,xn+1)*exx(yn+1,xn+1);
                
                
                % Interpolating squares of EPSxy=EPSyx from basic nodes
                % Horizontal EPSxy index
                xn=xnmin;
                % Vertical EPSxy index
                yn=ynmin;
                
                % Define and check normalized distances from marker to the upper left VX-node;
                dx=(MX(ym,xm)-gridx(xn))./xstp1(xn);
                dy=(MY(ym,xm)-gridy(yn))./ystp1(yn);
                
                % Calculate Marker EPSxy^2 from four surrounding nodes
                exym=0;
                exym=exym+(1.0-dx).*(1.0-dy).*exy(yn,xn)*exy(yn,xn);
                exym=exym+(1.0-dx).*dy.*exy(yn+1,xn)*exy(yn+1,xn);
                exym=exym+dx.*(1.0-dy).*exy(yn,xn+1)*exy(yn,xn+1);
                exym=exym+dx.*dy.*exy(yn+1,xn+1)*exy(yn+1,xn+1);
                %
                %Computing second strain rate invariant for the marker
                MEII(ym,xm)=(exxm+exym)^0.5;
            end
        end
    end


    % Moving Markers by velocity field
    if(markmove>0)
        % Create arrays for velocity of markers
        vxm=zeros(4,1);
        vym=zeros(4,1);
        % Marker cycle
        for xm = 1:1:mxnum
            for ym = 1:1:mynum  

                % Check markers inside the grid
                if (MX(ym,xm)>=0 && MX(ym,xm)<=xsize && MY(ym,xm)>=0 && MY(ym,xm)<=ysize) 

                    % Save marker coordinates
                    xcur=MX(ym,xm);
                    ycur=MY(ym,xm);
                    % Defining number of Runge-Kutta cycles
                    for rk=1:1:markmove

                        %  xn    V(xn,yn)--------------------V(xn+1,yn)
                        %           ?           ^                  ?
                        %           ?           ?                  ?
                        %           ?          dy                  ?
                        %           ?           ?                  ?
                        %           ?           v                  ?
                        %           ?<----dx--->o Mrho(xm,ym)       ?
                        %           ?                              ?
                        %           ?                              ?
                        %  xn+1  V(xn,yn+1)-------------------V(xn+1,yn+1)

                       % Define indexes for upper left BASIC node in the cell where the marker is
                        % using bisection
                        % Load horizontal and vertical indexes
                        if(rk==1)
                            xnmin=MXN(ym,xm);
                            ynmin=MYN(ym,xm);
                        else
                            % Find horizontal index
                            xnmin=1;
                            xnmax=xnum;
                            while ((xnmax-xnmin)>1)
                                % !!! SUBTRACT 0.5 since int16(0.5)=1
                                xn=double(int16((xnmax+xnmin)./2-0.5));
                                if(gridx(xn)>xcur)
                                    xnmax=xn;
                                else
                                    xnmin=xn;
                                end
                            end
                            % Check horizontal index
                            if (xnmin<1)
                                xnmin=1;
                            end
                            if (xnmin>xnum-1)
                                xnmin=xnum-1;
                            end
                            % Find vertical index
                            ynmin=1;
                            ynmax=ynum;
                            while ((ynmax-ynmin)>1)
                                % !!! SUBTRACT 0.5 since int16(0.5)=1
                                yn=double(int16((ynmax+ynmin)./2-0.5));
                                if(gridy(yn)>ycur)
                                    ynmax=yn;
                                else
                                    ynmin=yn;
                                end
                            end
                            % Check vertical index
                            if (ynmin<1)
                                ynmin=1;
                            end
                            if (ynmin>ynum-1)
                                ynmin=ynum-1;
                            end
                        end

                        % Define indexes for upper left node in the Vx-cell where the marker is
                        % Horizontal Vx index
                        xn=xnmin;
                        % Vertical Vx index
                        yn=ynmin;
                        if(ycur>gridcy(yn+1))
                            yn=yn+1;
                        end
                        if (yn>ynum)
                            yn=ynum;
                        end
                        
                        % Define and check normalized distances from marker to the upper left VX-node;
                        dx=(xcur-gridx(xn))./xstp1(xn);
                        dy=(ycur-gridcy(yn))./ystpc1(yn);
                        
                        % Calculate Marker velocity from four surrounding Vx nodes
                        vxm(rk)=0;
                        vxm(rk)=vxm(rk)+(1.0-dx).*(1.0-dy).*vx1(yn,xn);
                        vxm(rk)=vxm(rk)+(1.0-dx).*dy.*vx1(yn+1,xn);
                        vxm(rk)=vxm(rk)+dx.*(1.0-dy).*vx1(yn,xn+1);
                        vxm(rk)=vxm(rk)+dx.*dy.*vx1(yn+1,xn+1);

                        % Define indexes for upper left node in the VY-cell where the marker is
                        % Vertical Vy index
                        yn=ynmin;
                        % Horizontal Vy index
                        xn=xnmin;
                        if(xcur>gridcx(xn+1))
                            xn=xn+1;
                        end
                        if (xn>xnum)
                            xn=xnum;
                        end
                        
                        % Define and check normalized distances from marker to the upper left VX-node;
                        dx=(xcur-gridcx(xn))./xstpc1(xn);
                        dy=(ycur-gridy(yn))./ystp1(yn);
                        
                        % Calculate Marker velocity from four surrounding nodes
                        vym(rk)=0;
                        vym(rk)=vym(rk)+(1.0-dx).*(1.0-dy).*vy1(yn,xn);
                        vym(rk)=vym(rk)+(1.0-dx).*dy.*vy1(yn+1,xn);
                        vym(rk)=vym(rk)+dx.*(1.0-dy).*vy1(yn,xn+1);
                        vym(rk)=vym(rk)+dx.*dy.*vy1(yn+1,xn+1);

                        % Update coordinates for the next cycle
                        if(rk<4)
                            if (rk<3)
                                xcur=MX(ym,xm)+timestep/2*vxm(rk);
                                ycur=MY(ym,xm)+timestep/2*vym(rk);
                            else
                                xcur=MX(ym,xm)+timestep*vxm(rk);
                                ycur=MY(ym,xm)+timestep*vym(rk);
                            end
                        end
                    end
                    % Recompute velocity using 4-th order Runge_Kutta
                    if (markmove==4)
                        vxm(1)=(vxm(1)+2*vxm(2)+2*vxm(3)+vxm(4))/6;
                        vym(1)=(vym(1)+2*vym(2)+2*vym(3)+vym(4))/6;
                    end

                    % Displacing Marker according to its velocity
                    MX(ym,xm)=MX(ym,xm)+timestep*vxm(1);
                    MY(ym,xm)=MY(ym,xm)+timestep*vym(1);

    %                 % Recycling markers exiting the model
    %                 % Top->Bottom
    %                 if (MY(ym,xm)<0)
    %                     MY(ym,xm)=ysize+MY(ym,xm);
    %                     MTK(ym,xm)=MTK(ym,xm)-ttop+tbottom;
    %                 end
    %                  % Bottom->Top
    %                 if (MY(ym,xm)>ysize)
    %                     MY(ym,xm)=MY(ym,xm)-ysize;
    %                     MTK(ym,xm)=MTK(ym,xm)+ttop-tbottom;
    %                 end
    
                end
                
           end
        end
    end


    % Interpolating vx, vy for the basic grid
    vxb=zeros(ynum,xnum);
    vyb=zeros(ynum,xnum);
    for j=1:1:xnum
        for i=1:1:ynum
            vxb(i,j)=(vx1(i,j)+vx1(i+1,j))/2;
            vyb(i,j)=(vy1(i,j)+vy1(i,j+1))/2;
        end
    end
    
    figure(1);
    hold on;
    % Plot gridlines
    % Plot vertical gridlines
    for j=1:1:xnum
        plot([gridx(j)/1000 gridx(j)/1000],[0 ysize/1000],'k'); 
    end
    % Plot horizontal gridlines
    for j=1:1:ynum
        plot([0 xsize/1000],[gridy(j)/1000 gridy(j)/1000],'k'); 
    end
    hold off;
    axis ij image;
    box on;
%     print ('-dtiff', '-r300','-zbuffer ','fig1');


    
    figure(2);
    % Plotting isotherms
    contour(gridx/1000,gridy/1000,tk2,273:50:1273,'k');
    axis ij image;
    box on;
%     print ('-dtiff', '-r300','-zbuffer ','fig2');
    
    
    figure(3);
    % Plotting ToC
    pcolor(gridx/1000,gridy/1000,tk2-273);
%     colormap(Gray);
    shading interp;
    colorbar;
    title(['Step=',num2str(ntimestep),' Myr=',num2str(timesum*1e-6/(365.25*24*3600))]); 
    hold on;
    % Plotting velocity
    quiver(gridx(2:2:xnum-1)/1000,gridy(2:2:ynum-1)/1000,vxb(2:2:xnum-1,2:2:ynum-1),vyb(2:2:xnum-1,2:2:ynum-1),'k');
    hold off;
    axis ij image;
    box on;
%     print ('-dtiff', '-r300','-zbuffer ','fig3');

    
    figure(4);
    % Computing and plotting vrms, tavr
    vrms(ntimestep,1)=0;
    tavr(ntimestep,1)=0;
    timesum1(ntimestep)=timesum*1e-6/(365.25*24*3600);
    for j=1:1:xnum-1
        for i=1:1:ynum-1
            vrms(ntimestep,1)=vrms(ntimestep,1)+(vx1(i+1,j)^2+vx1(i+1,j+1)^2+vy1(i,j+1)^2+vy1(i+1,j+1)^2)/2*xstp1(j)*ystp1(i);
            tavr(ntimestep,1)=tavr(ntimestep,1)+(tk2(i,j)+tk2(i+1,j)+tk2(i,j+1)+tk2(i+1,j+1))/4*xstp1(j)*ystp1(i);
        end
    end
    vrms(ntimestep,1)=(vrms(ntimestep,1)/xsize/ysize)^0.5*ysize/(MKT(1)/MCP(1)/MRHO(1,1));
    tavr(ntimestep,1)=tavr(ntimestep,1)/xsize/ysize;
    subplot(3,1,1)
    plot(timesum1,tavr);
    title(['Tavr=',num2str(tavr(ntimestep,1))]);  

    subplot(3,1,2)
    plot(timesum1,vrms);
    title(['Vrms=',num2str(vrms(ntimestep,1))]);
    
    % Computing and plotting Nusselt number
    nnus1(ntimestep,1)=0;
    for j=1:1:xnum-1
        dtdy1=(tk2(1,j)-tk2(2,j)+tk2(1,j+1)-tk2(2,j+1))/2/ystp1(1);
        nnus1(ntimestep,1)=nnus1(ntimestep,1)-dtdy1*xstp1(j);
    end
    nnus1(ntimestep,1)=ysize*nnus1(ntimestep,1)/((tbottom-ttop)*xsize);
    subplot(3,1,3)
    plot(timesum1,nnus1,'k');
    title(['Nnus1=',num2str(nnus1(ntimestep,1))]);
%     print ('-dtiff', '-r300','-zbuffer ','fig4');


    
    timesum=timesum
    timestep=timestep
    
    % Advance in time
    timesum=timesum+timestep
    
    ntimestep=ntimestep

 

    pause(1)

end