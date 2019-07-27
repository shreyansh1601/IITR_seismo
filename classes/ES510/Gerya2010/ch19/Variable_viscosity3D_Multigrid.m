% Solution of Stokes and continuity equations in 3D with Multigrid
% based on V-cycle
% by using external function Stokes_smoother3D()
% Density distribution in the model corresponds to falling block test
% viscosity in the model is constant
% 
% Staggered Grid for Multigrid
% 
%     vx       vx       vx    
%
% vy  +---vy---+---vy---+   vy
%     |        |        |
%     vx   P   vx   P   vx    
%     |        |        |
% vy  +---vy---+---vy---+   vy
%     |        |        |
%     vx   P   vx   P   vx    
%     |        |        |
% vy  +---vy---+---vy---+   vy
%
%     vx       vx       vx    
% 
% Lines show basic grid
% Basic (density) nodes are shown with +
% Ghost nodes shown outside the basic grid
% are used for boundary conditions

% Clearing all variables and arrays
clear all;
% Clearing figures
clf;

% Model parameters
% Block density
rhoblock=3100.0;
% Block viscosity, Pa*s
etablock=1e+23;

% Medium density
rhomedium=3000.0;
% Medium viscosity, Pa*s
etamedium=1e+20;

% Acceleration of Gravity, m/s^2
g=9.81;
% Pressure in the upermost, lwftmost (first) cell
prfirst=0;
% Model size, m
xsize=100000.0;
ysize=100000.0;
zsize=100000.0;

% Multigrid parameters
% Number of resolution levels
levelnum=4;


% Defining limits for viscosity variations
etamin=etamedium;
etamax=etablock;
etamincur=etamedium;
etamaxcur=etamedium;
etaminstp=1;
etamaxstp=3.3333;
etacyc=15;
etannn=1;
resmax=1e-6;



% Iteration Parameters
% Total number of iteration cycles
inum=250;
% Relaxation coefficients for Gauss-Seidel iterations
% Stokes equations
krelaxs=0.9;
vkoef=1.0;
% Continuity equation
krelaxc=0.3;
pkoef=1.0;

% Number of smoothing iterations
iternum(1)=5;
iternum(2)=5*2;
iternum(3)=5*4;
iternum(4)=5*8;
iternum(5)=5*16;
iternum(6)=5*32;
iternum(7)=5*64;
iternum(8)=5*128;

% Defining resolutions on all levels
xnum(1)=25;
ynum(1)=25;
znum(1)=25;
xnum(2)=13;
ynum(2)=13;
znum(2)=13;
xnum(3)=7;
ynum(3)=7;
znum(3)=7;
xnum(4)=4;
ynum(4)=4;
znum(4)=4;

% xnum(1)=49;
% ynum(1)=49;
% znum(1)=49;
% xnum(2)=25;
% ynum(2)=25;
% znum(2)=25;
% xnum(3)=13;
% ynum(3)=13;
% znum(3)=13;
% xnum(4)=7;
% ynum(4)=7;
% znum(4)=7;


% xnum(1)=49;
% ynum(1)=49;
% znum(1)=49;
% xnum(2)=25;
% ynum(2)=25;
% znum(2)=25;
% xnum(3)=13;
% ynum(3)=13;
% znum(3)=13;
% xnum(4)=7;
% ynum(4)=7;
% znum(4)=7;
% xnum(5)=4;
% ynum(5)=4;
% znum(5)=4;

% xnum(1)=97;
% ynum(1)=97;
% znum(1)=97;
% xnum(2)=49;
% ynum(2)=49;
% znum(2)=49;
% xnum(3)=25;
% ynum(3)=25;
% znum(3)=25;
% xnum(4)=13;
% ynum(4)=13;
% znum(4)=13;
% xnum(5)=7;
% ynum(5)=7;
% znum(5)=7;
% xnum(6)=4;
% ynum(6)=4;
% znum(6)=4;


% Defining gridsteps on all levels
% Number 
for n=levelnum:-1:1
    % Grid steps
    xstp(n)=xsize./(xnum(n)-1);
    ystp(n)=ysize./(ynum(n)-1);
    zstp(n)=zsize./(znum(n)-1);
end


% FINEST (PRINCIPAL) GRID
% Defining density structure rho()
% Defining initial guesses for velocity vx() vy() and pressure pr()
% Computing right part of Stokes (RX, RY) and Continuity (RC) equation
% vx, vy, P
vx1=zeros(ynum(1)+1,xnum(1),znum(1)+1);
vy1=zeros(ynum(1),xnum(1)+1,znum(1)+1);
vz1=zeros(ynum(1)+1,xnum(1)+1,znum(1));
pr1=zeros(ynum(1)-1,xnum(1)-1,znum(1)-1);
% Right parts of equations
RX1=zeros(ynum(1)+1,xnum(1),znum(1)+1);
RY1=zeros(ynum(1),xnum(1)+1,znum(1)+1);
RZ1=zeros(ynum(1)+1,xnum(1)+1,znum(1));
RC1=zeros(ynum(1)-1,xnum(1)-1,znum(1)-1);
% Grid points cycle
for i=1:1:ynum(1);
    for j=1:1:xnum(1);
        for k=1:1:znum(1);
            
            % Density Structure
            % Relative distances for (i,j) basic node inside the grid
            dx=(j-1)*xstp(1)/xsize;
            dy=(i-1)*ystp(1)/ysize;
            dz=(k-1)*zstp(1)/zsize;
            % Density Structure (falling block)
            rho1(i,j,k)=rhomedium;
            if(dx>=0.4 && dx<=0.6 && dy>=0.4 && dy<=0.6 && dz>=0.4 && dz<=0.6)
                rho1(i,j,k)=rhoblock;
            end
            
            % Viscosity for normal stress (etan) for cells (in pressure nodes)
            if (i<ynum(1) && j<xnum(1) && k<znum(1))
                % Relative distances for (i,j) pressure nodes inside the grid
                dx=(j-0.5)*xstp(1)/xsize;
                dy=(i-0.5)*ystp(1)/ysize;
                dz=(k-0.5)*zstp(1)/zsize;
                etan1(i,j,k)=etamedium;
                if(dx>=0.4 && dx<=0.6 && dy>=0.4 && dy<=0.6 && dz>=0.4 && dz<=0.6)
                    etan1(i,j,k)=etablock;
                end
            end
            
            % Shear viscosity for SIGMAxy stress (etaxy) for z-directed edges
            if (k<znum(1))
                % Relative distances for (i,j) pressure nodes inside the grid
                dx=(j-1)*xstp(1)/xsize;
                dy=(i-1)*ystp(1)/ysize;
                dz=(k-0.5)*zstp(1)/zsize;
                etaxy1(i,j,k)=etamedium;
                if(dx>=0.4 && dx<=0.6 && dy>=0.4 && dy<=0.6 && dz>=0.4 && dz<=0.6)
                    etaxy1(i,j,k)=etablock;
                end
            end
            
            % Shear viscosity for SIGMAxz stress (etaxz) for y-directed edges
            if (i<ynum(1))
                % Relative distances for (i,j) pressure nodes inside the grid
                dx=(j-1)*xstp(1)/xsize;
                dy=(i-0.5)*ystp(1)/ysize;
                dz=(k-1)*zstp(1)/zsize;
                etaxz1(i,j,k)=etamedium;
                if(dx>=0.4 && dx<=0.6 && dy>=0.4 && dy<=0.6 && dz>=0.4 && dz<=0.6)
                    etaxz1(i,j,k)=etablock;
                end
            end
            
            % Shear viscosity for SIGMAyz stress (etayz) for x-directed edges
            if (j<xnum(1))
                % Relative distances for (i,j) pressure nodes inside the grid
                dx=(j-0.5)*xstp(1)/xsize;
                dy=(i-1)*ystp(1)/ysize;
                dz=(k-1)*zstp(1)/zsize;
                etayz1(i,j,k)=etamedium;
                if(dx>=0.4 && dx<=0.6 && dy>=0.4 && dy<=0.6 && dz>=0.4 && dz<=0.6)
                    etayz1(i,j,k)=etablock;
                end
            end

            % Pressure update in internal nodes (lithostatic pressure)
            if(i<ynum(1) && j>1 && k>1)
                if(i==1)
                    pr1(i,j-1,k-1)=prfirst;
                else
                    pr1(i,j-1,k-1)=pr1(i-1,j-1,k-1)+ystp(1)*g*(rho1(i,j,k)+rho1(i,j-1,k)+rho1(i,j,k-1)+rho1(i,j-1,k-1))/4;
                end
            end
            
            % Right part of x-Stokes Equation
            if(j>1 && i>1 && k>1 && j<xnum(1))
                RX1(i,j,k)=0;
            end
            
            % Right part of y-Stokes Equation
            if(j>1 && k>1 && i>1 && i<ynum(1))
                RY1(i,j,k)=-g*(rho1(i,j,k)+rho1(i,j-1,k)+rho1(i,j,k-1)+rho1(i,j-1,k-1))/4;
            end
            
            % Right part of z-Stokes Equation
            if(j>1 && i>1 && k>1 && k<znum(1))
                RZ1(i,j,k)=0;
            end
            
        end
    end
end


% Interpolation (restriction) of viscosity to coarser levels
for n=1:1:levelnum-1;
    % Action depends on the level of resolution
    switch n
            
        % Level 1 (principal grid)
        case 1
        if (levelnum>n)
            [etaxy2 etaxz2 etayz2 etan2]=Viscosity_restriction3D(n,xnum,ynum,znum,xstp,ystp,zstp,etaxy1,etaxz1,etayz1,etan1);
        end
        % Level 2 
        case 2
        if (levelnum>n)
            [etaxy3 etaxz3 etayz3 etan3]=Viscosity_restriction3D(n,xnum,ynum,znum,xstp,ystp,zstp,etaxy2,etaxz2,etayz2,etan2);
        end
        % Level 3 
        case 3
        if (levelnum>n)
            [etaxy4 etaxz4 etayz4 etan4]=Viscosity_restriction3D(n,xnum,ynum,znum,xstp,ystp,zstp,etaxy3,etaxz3,etayz3,etan3);
        end
        % Level 4 
        case 4
        if (levelnum>n)
            [etaxy5 etaxz5 etayz5 etan5]=Viscosity_restriction3D(n,xnum,ynum,znum,xstp,ystp,zstp,etaxy4,etaxz4,etayz4,etan4);
        end
        % Level 5 
        case 5
        if (levelnum>n)
            [etaxy6 etaxz6 etayz6 etan6]=Viscosity_restriction3D(n,xnum,ynum,znum,xstp,ystp,zstp,etaxy5,etaxz5,etayz5,etan5);
        end
        % Level 6 
        case 6
        if (levelnum>n)
            [etaxy7 etaxz7 etayz7 etan7]=Viscosity_restriction3D(n,xnum,ynum,znum,xstp,ystp,zstp,etaxy6,etaxz6,etayz6,etan6);
        end
        % Level 7 
        case 7
        if (levelnum>n)
            [etaxy8 etaxz8 etayz8 etan8]=Viscosity_restriction3D(n,xnum,ynum,znum,xstp,ystp,zstp,etaxy7,etaxz7,etayz7,etan7);
        end
        % No further restriction from last Level 8
    end
end




% Saving viscosity
if (levelnum>0) 
    etaxyo1=etaxy1; 
    etaxzo1=etaxz1; 
    etayzo1=etayz1; 
    etano1=etan1; 
end
if (levelnum>1) 
    etaxyo2=etaxy2; 
    etaxzo2=etaxz2; 
    etayzo2=etayz2; 
    etano2=etan2; 
end
if (levelnum>2) 
    etaxyo3=etaxy3; 
    etaxzo3=etaxz3; 
    etayzo3=etayz3; 
    etano3=etan3; 
end
if (levelnum>3) 
    etaxyo4=etaxy4; 
    etaxzo4=etaxz4; 
    etayzo4=etayz4; 
    etano4=etan4; 
end
if (levelnum>4) 
    etaxyo5=etaxy5; 
    etaxzo5=etaxz5; 
    etayzo5=etayz5; 
    etano5=etan5; 
end
if (levelnum>5) 
    etaxyo6=etaxy6; 
    etaxzo6=etaxz6; 
    etayzo6=etayz6; 
    etano6=etan6; 
end
if (levelnum>6) 
    etaxyo7=etaxy7; 
    etaxzo7=etaxz7; 
    etayzo7=etayz7; 
    etano7=etan7; 
end
if (levelnum>7) 
    etaxyo8=etaxy8; 
    etaxzo8=etaxz8; 
    etayzo8=etayz8; 
    etano8=etan8; 
end


% Renormalizing viscosity
if (etamaxcur<etamax || etamincur>etamin)
    % Normalization coefficient
    etakf=log(etamaxcur/etamincur)/log(etamax/etamin);
    % Recomputing viscosity
    if (levelnum>0) 
        etaxy1=etamincur*exp(log(etaxyo1/etamin)*etakf);
        etaxz1=etamincur*exp(log(etaxzo1/etamin)*etakf);
        etayz1=etamincur*exp(log(etayzo1/etamin)*etakf);
        etan1=etamincur*exp(log(etano1/etamin)*etakf);
    end
    if (levelnum>1) 
        etaxy2=etamincur*exp(log(etaxyo2/etamin)*etakf);
        etaxz2=etamincur*exp(log(etaxzo2/etamin)*etakf);
        etayz2=etamincur*exp(log(etayzo2/etamin)*etakf);
        etan2=etamincur*exp(log(etano2/etamin)*etakf);
    end
    if (levelnum>2) 
        etaxy3=etamincur*exp(log(etaxyo3/etamin)*etakf);
        etaxz3=etamincur*exp(log(etaxzo3/etamin)*etakf);
        etayz3=etamincur*exp(log(etayzo3/etamin)*etakf);
        etan3=etamincur*exp(log(etano3/etamin)*etakf);
    end
    if (levelnum>3) 
        etaxy4=etamincur*exp(log(etaxyo4/etamin)*etakf);
        etaxz4=etamincur*exp(log(etaxzo4/etamin)*etakf);
        etayz4=etamincur*exp(log(etayzo4/etamin)*etakf);
        etan4=etamincur*exp(log(etano4/etamin)*etakf);
    end
    if (levelnum>4) 
        etaxy5=etamincur*exp(log(etaxyo5/etamin)*etakf);
        etaxz5=etamincur*exp(log(etaxzo5/etamin)*etakf);
        etayz5=etamincur*exp(log(etayzo5/etamin)*etakf);
        etan5=etamincur*exp(log(etano5/etamin)*etakf);
    end
    if (levelnum>5) 
        etaxy6=etamincur*exp(log(etaxyo6/etamin)*etakf);
        etaxz6=etamincur*exp(log(etaxzo6/etamin)*etakf);
        etayz6=etamincur*exp(log(etayzo6/etamin)*etakf);
        etan6=etamincur*exp(log(etano6/etamin)*etakf);
    end
    if (levelnum>6) 
        etaxy7=etamincur*exp(log(etaxyo7/etamin)*etakf);
        etaxz7=etamincur*exp(log(etaxzo7/etamin)*etakf);
        etayz7=etamincur*exp(log(etayzo7/etamin)*etakf);
        etan7=etamincur*exp(log(etano7/etamin)*etakf);
    end
    if (levelnum>7) 
        etaxy8=etamincur*exp(log(etaxyo8/etamin)*etakf);
        etaxz8=etamincur*exp(log(etaxzo8/etamin)*etakf);
        etayz8=etamincur*exp(log(etayzo8/etamin)*etakf);
        etan8=etamincur*exp(log(etano8/etamin)*etakf);
    end
end




% Main Multigrid cycle
% We put 8 levels (for the case if we want to increase resolution)
% We use different arrays for different levels 
% since grid resolutions is different
% Multigrid schedule = V-cycle 
for niter=1:1:inum;
    
    
    % Smoothing+restriction cycle
    for n=1:1:levelnum;
        % Action depends on the level of resolution
        switch n
            
            % Leve 1 (principal grid)
            case 1
            % Smoothing operation: solving of Stokes and Continuity equations on nodes
            % and computing residuals
            % by calling function Stokes_smoother()
            [vx1,resx1,vy1,resy1,vz1,resz1,pr1,resc1]=Stokes_Continuity3D_viscous_smoother(prfirst,etaxy1,etaxz1,etayz1,etan1,iternum(n),krelaxs,krelaxc,xnum(n),ynum(n),znum(n),xstp(n),ystp(n),zstp(n),RX1,vx1,RY1,vy1,RZ1,vz1,RC1,pr1);
            % Restriction operation: 
            % Interpolating residuals to coarcer level (k+1) 
            % to produce right parts for this coarcer level
            if (levelnum>n)
                [RX2,RY2,RZ2,RC2]=Stokes_Continuity3D_viscous_restriction(n,xnum,ynum,znum,xstp,ystp,zstp,resx1,resy1,resz1,resc1,etan1,etan2);
            end
        
            % Level 2
            case 2
            %Making initial approximation for vx,vy,vz, P corrections (zeros)
            vx2=zeros(ynum(n)+1,xnum(n),znum(n)+1);
            vy2=zeros(ynum(n),xnum(n)+1,znum(n)+1);
            vz2=zeros(ynum(n)+1,xnum(n)+1,znum(n));
            pr2=zeros(ynum(n)-1,xnum(n)-1,znum(n)-1);
            % Smoothing operation: 
            [vx2,resx2,vy2,resy2,vz2,resz2,pr2,resc2]=Stokes_Continuity3D_viscous_smoother(0,etaxy2,etaxz2,etayz2,etan2,iternum(n),krelaxs,krelaxc,xnum(n),ynum(n),znum(n),xstp(n),ystp(n),zstp(n),RX2,vx2,RY2,vy2,RZ2,vz2,RC2,pr2);
            % Restriction operation: 
            if (levelnum>n)
                [RX3,RY3,RZ3,RC3]=Stokes_Continuity3D_viscous_restriction(n,xnum,ynum,znum,xstp,ystp,zstp,resx2,resy2,resz2,resc2,etan2,etan3);
            end
        
            % Level 3
            case 3
            %Making initial approximation for vx,vy,vz, P corrections (zeros)
            vx3=zeros(ynum(n)+1,xnum(n),znum(n)+1);
            vy3=zeros(ynum(n),xnum(n)+1,znum(n)+1);
            vz3=zeros(ynum(n)+1,xnum(n)+1,znum(n));
            pr3=zeros(ynum(n)-1,xnum(n)-1,znum(n)-1);
            % Smoothing operation: 
            [vx3,resx3,vy3,resy3,vz3,resz3,pr3,resc3]=Stokes_Continuity3D_viscous_smoother(0,etaxy3,etaxz3,etayz3,etan3,iternum(n),krelaxs,krelaxc,xnum(n),ynum(n),znum(n),xstp(n),ystp(n),zstp(n),RX3,vx3,RY3,vy3,RZ3,vz3,RC3,pr3);
            % Restriction operation: 
            if (levelnum>n)
                [RX4,RY4,RZ4,RC4]=Stokes_Continuity3D_viscous_restriction(n,xnum,ynum,znum,xstp,ystp,zstp,resx3,resy3,resz3,resc3,etan3,etan4);
            end
            
            % Level 4
            case 4
            %Making initial approximation for vx,vy,vz, P corrections (zeros)
            vx4=zeros(ynum(n)+1,xnum(n),znum(n)+1);
            vy4=zeros(ynum(n),xnum(n)+1,znum(n)+1);
            vz4=zeros(ynum(n)+1,xnum(n)+1,znum(n));
            pr4=zeros(ynum(n)-1,xnum(n)-1,znum(n)-1);
            % Smoothing operation: 
            [vx4,resx4,vy4,resy4,vz4,resz4,pr4,resc4]=Stokes_Continuity3D_viscous_smoother(0,etaxy4,etaxz4,etayz4,etan4,iternum(n),krelaxs,krelaxc,xnum(n),ynum(n),znum(n),xstp(n),ystp(n),zstp(n),RX4,vx4,RY4,vy4,RZ4,vz4,RC4,pr4);
            % Restriction operation: 
            if (levelnum>n)
                [RX5,RY5,RZ5,RC5]=Stokes_Continuity3D_viscous_restriction(n,xnum,ynum,znum,xstp,ystp,zstp,resx4,resy4,resz4,resc4,etan4,etan5);
            end
            
            % Level 5
            case 5
            %Making initial approximation for vx,vy,vz, P corrections (zeros)
            vx5=zeros(ynum(n)+1,xnum(n),znum(n)+1);
            vy5=zeros(ynum(n),xnum(n)+1,znum(n)+1);
            vz5=zeros(ynum(n)+1,xnum(n)+1,znum(n));
            pr5=zeros(ynum(n)-1,xnum(n)-1,znum(n)-1);
            % Smoothing operation: 
            [vx5,resx5,vy5,resy5,vz5,resz5,pr5,resc5]=Stokes_Continuity3D_viscous_smoother(0,etaxy5,etaxz5,etayz5,etan5,iternum(n),krelaxs,krelaxc,xnum(n),ynum(n),znum(n),xstp(n),ystp(n),zstp(n),RX5,vx5,RY5,vy5,RZ5,vz5,RC5,pr5);
            % Restriction operation: 
            if (levelnum>n)
                [RX6,RY6,RZ6,RC6]=Stokes_Continuity3D_viscous_restriction(n,xnum,ynum,znum,xstp,ystp,zstp,resx5,resy5,resz5,resc5,etan5,etan6);
            end
             
            % Level 6
            case 6
            %Making initial approximation for vx,vy,vz, P corrections (zeros)
            vx6=zeros(ynum(n)+1,xnum(n),znum(n)+1);
            vy6=zeros(ynum(n),xnum(n)+1,znum(n)+1);
            vz6=zeros(ynum(n)+1,xnum(n)+1,znum(n));
            pr6=zeros(ynum(n)-1,xnum(n)-1,znum(n)-1);
            % Smoothing operation: 
            [vx6,resx6,vy6,resy6,vz6,resz6,pr6,resc6]=Stokes_Continuity3D_viscous_smoother(0,etaxy6,etaxz6,etayz6,etan6,iternum(n),krelaxs,krelaxc,xnum(n),ynum(n),znum(n),xstp(n),ystp(n),zstp(n),RX6,vx6,RY6,vy6,RZ6,vz6,RC6,pr6);
            % Restriction operation: 
            if (levelnum>n)
                [RX7,RY7,RZ7,RC7]=Stokes_Continuity3D_viscous_restriction(n,xnum,ynum,znum,xstp,ystp,zstp,resx6,resy6,resz6,resc6,etan6,etan7);
            end
             
            % Level 7
            case 7
            %Making initial approximation for vx,vy,vz, P corrections (zeros)
            vx7=zeros(ynum(n)+1,xnum(n),znum(n)+1);
            vy7=zeros(ynum(n),xnum(n)+1,znum(n)+1);
            vz7=zeros(ynum(n)+1,xnum(n)+1,znum(n));
            pr7=zeros(ynum(n)-1,xnum(n)-1,znum(n)-1);
            % Smoothing operation: 
            [vx7,resx7,vy7,resy7,vz7,resz7,pr7,resc7]=Stokes_Continuity3D_viscous_smoother(0,etaxy7,etaxz7,etayz7,etan7,iternum(n),krelaxs,krelaxc,xnum(n),ynum(n),znum(n),xstp(n),ystp(n),zstp(n),RX7,vx7,RY7,vy7,RZ7,vz7,RC7,pr7);
            % Restriction operation: 
            if (levelnum>n)
                [RX8,RY8,RZ8,RC8]=Stokes_Continuity3D_viscous_restriction(n,xnum,ynum,znum,xstp,ystp,zstp,resx7,resy7,resz7,resc7,etan7,etan8);
            end
             
            % Level 8
            case 8
            %Making initial approximation for vx,vy,vz, P corrections (zeros)
            vx8=zeros(ynum(n)+1,xnum(n),znum(n)+1);
            vy8=zeros(ynum(n),xnum(n)+1,znum(n)+1);
            vz8=zeros(ynum(n)+1,xnum(n)+1,znum(n));
            pr8=zeros(ynum(n)-1,xnum(n)-1,znum(n)-1);
            % Smoothing operation: 
            [vx8,resx8,vy8,resy8,vz8,resz8,pr8,resc8]=Stokes_Continuity3D_viscous_smoother(0,etaxy8,etaxz8,etayz8,etan8,iternum(n),krelaxs,krelaxc,xnum(n),ynum(n),znum(n),xstp(n),ystp(n),zstp(n),RX8,vx8,RY8,vy8,RZ8,vz8,RC8,pr8);
            % No Restriction operation on this last level 

        end
    end
   
    


    
    
    % Smoothing+prolongation cycle
    for n=levelnum:-1:1;
        % Action depends on the level of resolution
        switch n

            % Leve 1 (principal grid)
            case 1
            % Smoothing operation: solving of Stokes and Continuity equations on nodes
            % and computing residuals
            % by calling function Stokes_smoother()
            [vx1,resx1,vy1,resy1,vz1,resz1,pr1,resc1]=Stokes_Continuity3D_viscous_smoother(prfirst,etaxy1,etaxz1,etayz1,etan1,iternum(n),krelaxs,krelaxc,xnum(n),ynum(n),znum(n),xstp(n),ystp(n),zstp(n),RX1,vx1,RY1,vy1,RZ1,vz1,RC1,pr1);
            % No prolongation operation on finest level
        
            % Level 2
            case 2
            % Smoothing operation: 
            [vx2,resx2,vy2,resy2,vz2,resz2,pr2,resc2]=Stokes_Continuity3D_viscous_smoother(0,etaxy2,etaxz2,etayz2,etan2,iternum(n),krelaxs,krelaxc,xnum(n),ynum(n),znum(n),xstp(n),ystp(n),zstp(n),RX2,vx2,RY2,vy2,RZ2,vz2,RC2,pr2);
            % Prolongation operation: 
            % Interpolating dvx,dvy,dpr corrections to finer level (k+1) 
            % and computing solution update for this finer level
            [dvx1,dvy1,dvz1,dpr1]=Stokes_Continuity3D_prolongation(n,xnum,ynum,znum,xstp,ystp,zstp,vx2,vy2,vz2,pr2);
            vx1=vx1+dvx1*vkoef;
            vy1=vy1+dvy1*vkoef;
            vz1=vz1+dvz1*vkoef;
            pr1=pr1+dpr1*pkoef;
        
            % Level 3
            case 3
            % Smoothing operation: 
            [vx3,resx3,vy3,resy3,vz3,resz3,pr3,resc3]=Stokes_Continuity3D_viscous_smoother(0,etaxy3,etaxz3,etayz3,etan3,iternum(n),krelaxs,krelaxc,xnum(n),ynum(n),znum(n),xstp(n),ystp(n),zstp(n),RX3,vx3,RY3,vy3,RZ3,vz3,RC3,pr3);
            % Prolongation operation: 
            [dvx2,dvy2,dvz2,dpr2]=Stokes_Continuity3D_prolongation(n,xnum,ynum,znum,xstp,ystp,zstp,vx3,vy3,vz3,pr3);
            vx2=vx2+dvx2*vkoef;
            vy2=vy2+dvy2*vkoef;
            vz2=vz2+dvz2*vkoef;
            pr2=pr2+dpr2*pkoef;
            
            % Level 4
            case 4
            % Smoothing operation: 
            [vx4,resx4,vy4,resy4,vz4,resz4,pr4,resc4]=Stokes_Continuity3D_viscous_smoother(0,etaxy4,etaxz4,etayz4,etan4,iternum(n),krelaxs,krelaxc,xnum(n),ynum(n),znum(n),xstp(n),ystp(n),zstp(n),RX4,vx4,RY4,vy4,RZ4,vz4,RC4,pr4);
            % Prolongation operation: 
            [dvx3,dvy3,dvz3,dpr3]=Stokes_Continuity3D_prolongation(n,xnum,ynum,znum,xstp,ystp,zstp,vx4,vy4,vz4,pr4);
            vx3=vx3+dvx3*vkoef;
            vy3=vy3+dvy3*vkoef;
            vz3=vz3+dvz3*vkoef;
            pr3=pr3+dpr3*pkoef;
            
             % Level 5
            case 5
            % Smoothing operation: 
            [vx5,resx5,vy5,resy5,vz5,resz5,pr5,resc5]=Stokes_Continuity3D_viscous_smoother(0,etaxy5,etaxz5,etayz5,etan5,iternum(n),krelaxs,krelaxc,xnum(n),ynum(n),znum(n),xstp(n),ystp(n),zstp(n),RX5,vx5,RY5,vy5,RZ5,vz5,RC5,pr5);
            % Prolongation operation: 
            [dvx4,dvy4,dvz4,dpr4]=Stokes_Continuity3D_prolongation(n,xnum,ynum,znum,xstp,ystp,zstp,vx5,vy5,vz5,pr5);
            vx4=vx4+dvx4*vkoef;
            vy4=vy4+dvy4*vkoef;
            vz4=vz4+dvz4*vkoef;
            pr4=pr4+dpr4*pkoef;
            
            
            % Level 6
            case 6
            % Smoothing operation: 
            [vx6,resx6,vy6,resy6,vz6,resz6,pr6,resc6]=Stokes_Continuity3D_viscous_smoother(0,etaxy6,etaxz6,etayz6,etan6,iternum(n),krelaxs,krelaxc,xnum(n),ynum(n),znum(n),xstp(n),ystp(n),zstp(n),RX6,vx6,RY6,vy6,RZ6,vz6,RC6,pr6);
            % Prolongation operation: 
            [dvx5,dvy5,dvz5,dpr5]=Stokes_Continuity3D_prolongation(n,xnum,ynum,znum,xstp,ystp,zstp,vx6,vy6,vz6,pr6);
            vx5=vx5+dvx5*vkoef;
            vy5=vy5+dvy5*vkoef;
            vz5=vz5+dvz5*vkoef;
            pr5=pr5+dpr5*pkoef;
            
            % Level 7
            case 7
            % Smoothing operation: 
            [vx7,resx7,vy7,resy7,vz7,resz7,pr7,resc7]=Stokes_Continuity3D_viscous_smoother(0,etaxy7,etaxz7,etayz7,etan7,iternum(n),krelaxs,krelaxc,xnum(n),ynum(n),znum(n),xstp(n),ystp(n),zstp(n),RX7,vx7,RY7,vy7,RZ7,vz7,RC7,pr7);
            % Prolongation operation: 
            [dvx6,dvy6,dvz6,dpr6]=Stokes_Continuity3D_prolongation(n,xnum,ynum,znum,xstp,ystp,zstp,vx7,vy7,vz7,pr7);
            vx6=vx6+dvx6*vkoef;
            vy6=vy6+dvy6*vkoef;
            vz6=vz6+dvz6*vkoef;
            pr6=pr6+dpr6*pkoef;
            
            % Level 8
            case 8
            % Smoothing operation: 
            [vx8,resx8,vy8,resy8,vz8,resz8,pr8,resc8]=Stokes_Continuity3D_viscous_smoother(0,etaxy8,etaxz8,etayz8,etan8,iternum(n),krelaxs,krelaxc,xnum(n),ynum(n),znum(n),xstp(n),ystp(n),zstp(n),RX8,vx8,RY8,vy8,RZ8,vz8,RC8,pr8);
            % Prolongation operation: 
            [dvx7,dvy7,dvz7,dpr7]=Stokes_Continuity3D_prolongation(n,xnum,ynum,znum,xstp,ystp,zstp,vx8,vy8,vz8,pr8);
            vx7=vx7+dvx7*vkoef;
            vy7=vy7+dvy7*vkoef;
            vz7=vz7+dvz7*vkoef;
            pr7=pr7+dpr7*pkoef;
        end
    end 
 
    % Defining scale for Stokes residuals from y-Stokes equation
    % dSIGMAij/dj-dP/di=-RHO*gi=0  => Stokes scale=abs(RHO*gi)
    stokesscale= rhoblock*g;
    % Defining scale for Continuity residuals from y-Stokes equation
    % dvx/dx+dvy/dy=0 can be transformed to 2ETA(dvx/dx+dvy/dy)/dx=0 
    % which is similar to dSIGMAij/dj and has scale given above
    % therefore continuity scale = scale=abs(RHO*gi/ETA*dx)
    continscale= rhoblock*g/etamedium*xstp(1);


    % Plotting Residuals for x-Stokes as surface
    figure(1);
    subplot(2,3,1)
    resx0=resx1/stokesscale;
    surf(resx0(:,:,(znum(1)-1)/2+1));
    shading interp;
    light;
    lighting phong;
    axis tight;
    zlabel('residual x-Stokes')
%     title(['x-Stokes residual, Gauss-Seidel cycle = ',num2str(niter)])
%     colorbar;
    
    % Plotting Residuals for y-Stokes as surface
    subplot(2,3,2)
    resy0=resy1/stokesscale;
    surf(resy0(:,:,(znum(1)-1)/2+1));
    shading interp;
    light;
    lighting phong;
    axis tight;
    zlabel('residual Y-stokes')
%     title(['y-Stokes residual, Gauss-Seidel cycle = ',num2str(niter)])
%     colorbar;

    % Plotting Residuals for z-Stokes as surface
    subplot(2,3,3)
    resz0=resz1/stokesscale;
    surf(resz0(:,:,(znum(1)-1)/2+1));
    %surf(resz2(:,:,(znum(1)-1)/4+1));
    shading interp;
    light;
    lighting phong;
    axis tight;
    zlabel('residual Z-stokes')
%     title(['y-Stokes residual, Gauss-Seidel cycle = ',num2str(niter)])
%     colorbar;

    % Plotting Residuals for Continuity as surface
    subplot(2,3,4)
     resc0=resc1/continscale;
    surf(resc0(:,:,(znum(1)-1)/2+1));
    shading interp;
    light;
    lighting phong;
    axis tight;
    zlabel('residual continuity')
%     title(['Continuity residual, Gauss-Seidel cycle = ',num2str(niter)])
%     colorbar;
    
    % Computing mean square residuals
    resx00(niter)=0;
    resy00(niter)=0;
    resz00(niter)=0;
    resc00(niter)=0;
    for i=1:1:ynum(1)
        for j=1:1:xnum(1)
            for k=1:1:znum(1)
                % x-Stokes
                if (i>1 && j>1 && k>1 && j<xnum(1))
                    resx00(niter)=resx00(niter)+resx0(i,j,k)^2;
                end
                % y-Stokes
                if (i>1 && j>1 && k>1 && i<ynum(1))
                    resy00(niter)=resy00(niter)+resy0(i,j,k)^2;
                end
                % z-Stokes
                if (i>1 && j>1 && k>1 && k<znum(1))
                    resz00(niter)=resz00(niter)+resz0(i,j,k)^2;
                end
                % Continuity
                if (j<xnum(1) && i<ynum(1) && k<znum(1))
                    resc00(niter)=resc00(niter)+resc0(i,j,k)^2;
                end
            end
        end
    end
    resx00(niter)=log10((resx00(niter)/((ynum(1)-1)*(xnum(1)-2)*(znum(1)-1)))^0.5);
    resy00(niter)=log10((resy00(niter)/((ynum(1)-2)*(xnum(1)-1)*(znum(1)-1)))^0.5);
    resz00(niter)=log10((resz00(niter)/((ynum(1)-1)*(xnum(1)-1)*(znum(1)-2)))^0.5);
    resc00(niter)=log10((resc00(niter)/((ynum(1)-1)*(xnum(1)-1)*(znum(1)-1)))^0.5);
   

    % Plotting Mean residuals
    subplot(2,3,5)
    plot(resc00, 'k');
    hold on
    plot(resx00, 'b');
    hold on
    plot(resy00, 'g');
    hold on
    plot(resz00, 'r');
    xlabel('V-cycles')
    ylabel('log(Residuals)')

    zlabel(['RMS errors, cycle ',num2str(niter)])
%     title(['RMS errors, Multigrid V-cycle = ',num2str(niter)])
%     colorbar;

 
    % Plotting Residuals for Continuity as surface
    subplot(2,3,6)
    surf(vy1(:,:,(znum(1)-1)/2+1));
    shading interp;
    light;
    lighting phong;
    axis tight;
    zlabel('Vy at z=0.5')
%     title(['Continuity residual, Gauss-Seidel cycle = ',num2str(niter)])
%     colorbar;

    pause(.1)

    % Renormalizing viscosity
    % Add cycle counter
    etannn=etannn+1;
    restot=(10^resx00(niter)+10^resx00(niter)+10^resx00(niter)+10^resx00(niter));
    if((etannn>etacyc || restot<resmax) && (etamaxcur<etamax || etamincur>etamin))
        % Reset counter
        etannn=1;
        % Change lower viscosity limit
        etamincur=etamincur*etaminstp;
        if (etamincur<etamin)
            etamincur=etamin; 
        end
        % Change upper viscosity limit
        etamaxcur=etamaxcur*etamaxstp;
        if (etamaxcur>etamax)
            etamaxcur=etamax; 
        end
        % Normalization coefficient
        etakf=log(etamaxcur/etamincur)/log(etamax/etamin);
        % Recomputing viscosity
        if (levelnum>0) 
            etaxy1=etamincur*exp(log(etaxyo1/etamin)*etakf);
            etaxz1=etamincur*exp(log(etaxzo1/etamin)*etakf);
            etayz1=etamincur*exp(log(etayzo1/etamin)*etakf);
            etan1=etamincur*exp(log(etano1/etamin)*etakf);
        end
        if (levelnum>1) 
            etaxy2=etamincur*exp(log(etaxyo2/etamin)*etakf);
            etaxz2=etamincur*exp(log(etaxzo2/etamin)*etakf);
            etayz2=etamincur*exp(log(etayzo2/etamin)*etakf);
            etan2=etamincur*exp(log(etano2/etamin)*etakf);
        end
        if (levelnum>2) 
            etaxy3=etamincur*exp(log(etaxyo3/etamin)*etakf);
            etaxz3=etamincur*exp(log(etaxzo3/etamin)*etakf);
            etayz3=etamincur*exp(log(etayzo3/etamin)*etakf);
            etan3=etamincur*exp(log(etano3/etamin)*etakf);
        end
        if (levelnum>3) 
            etaxy4=etamincur*exp(log(etaxyo4/etamin)*etakf);
            etaxz4=etamincur*exp(log(etaxzo4/etamin)*etakf);
            etayz4=etamincur*exp(log(etayzo4/etamin)*etakf);
            etan4=etamincur*exp(log(etano4/etamin)*etakf);
        end
        if (levelnum>4) 
            etaxy5=etamincur*exp(log(etaxyo5/etamin)*etakf);
            etaxz5=etamincur*exp(log(etaxzo5/etamin)*etakf);
            etayz5=etamincur*exp(log(etayzo5/etamin)*etakf);
            etan5=etamincur*exp(log(etano5/etamin)*etakf);
        end
        if (levelnum>5) 
            etaxy6=etamincur*exp(log(etaxyo6/etamin)*etakf);
            etaxz6=etamincur*exp(log(etaxzo6/etamin)*etakf);
            etayz6=etamincur*exp(log(etayzo6/etamin)*etakf);
            etan6=etamincur*exp(log(etano6/etamin)*etakf);
        end
        if (levelnum>6) 
            etaxy7=etamincur*exp(log(etaxyo7/etamin)*etakf);
            etaxz7=etamincur*exp(log(etaxzo7/etamin)*etakf);
            etayz7=etamincur*exp(log(etayzo7/etamin)*etakf);
            etan7=etamincur*exp(log(etano7/etamin)*etakf);
        end
        if (levelnum>7) 
            etaxy8=etamincur*exp(log(etaxyo8/etamin)*etakf);
            etaxz8=etamincur*exp(log(etaxzo8/etamin)*etakf);
            etayz8=etamincur*exp(log(etayzo8/etamin)*etakf);
            etan8=etamincur*exp(log(etano8/etamin)*etakf);
        end
    end


    
end

