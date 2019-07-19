% Solution of Stokes and continuity equations 
% with variable viscosity in 2D with Multigrid
% based on V-cycle
% by using external function Stokes_Continuity_viscous_smoother()
% Density distribution in the model corresponds to falling block test
% viscosity in the model is variable
% Setup corresponds to falling block test
% Grid resolutions for different levels are defined independently
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

% % Model parameters
% Block density
rhoblock=3100.0;
% Block viscosity, Pa*s
etablock=1e+26;

% Medium density
rhomedium=3000.0;
% Medium viscosity, Pa*s
etamedium=1e+20;

% Defining limits for viscosity variations
etamin=etamedium;
etamax=etablock;
etaminbeg=etamedium;
etamaxbeg=etamedium;
etamaxstp=10.0001;
etacyc=3;
etannn=1;
etamincur=etaminbeg;
etamaxcur=etamaxbeg;

% Acceleration of Gravity, m/s^2
g=9.81;
% Pressure in the upermost, lwftmost (first) cell
prfirst=0;
% Model size, m
xsize=100000.0;
ysize=100000.0;

% Multigrid parameters
% Way of avearing viscosity for interpolating to coarcer levels
% 0 - arithmetic
% 1 - geometric
% 2 - harmonic
etamean=0;
% Number of resolution levels
levelnum=5;
% Iteration Parameters
% Total number of iteration cycles
inum=500;
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


% Defining resolutions on different levels
xnum(1)=49;
ynum(1)=49;
xnum(2)=31;
ynum(2)=31;
xnum(3)=21;
ynum(3)=21;
xnum(4)=13;
ynum(4)=13;
xnum(5)=7;
ynum(5)=7;
xnum(6)=7;
ynum(6)=7;
xnum(7)=7;
ynum(7)=7;
xnum(8)=7;
ynum(8)=7;


% Defining gridsteps on all levels
% Number 
for k=levelnum:-1:1
    % Grid steps
    xstp(k)=xsize./(xnum(k)-1);
    ystp(k)=ysize./(ynum(k)-1);
end

% FINEST (PRINCIPAL) GRID
% Defining density rho() viscosity for shear stress (etas) and 
% viscosity for normal stress (etan) (in cells)
% Defining initial guesses for velocity vx() vy() and pressure pr()
% Computing right part of Stokes (RX, RY) and Continuity (RC) equation
% vx, vy, P
vx1=zeros(ynum(1)+1,xnum(1));
vy1=zeros(ynum(1),xnum(1)+1);
pr1=zeros(ynum(1)-1,xnum(1)-1);
vx0=vx1;
vy0=vy1;
pr0=pr1;
% Right parts of equations
RX0=zeros(ynum(1)+1,xnum(1));
RY0=zeros(ynum(1),xnum(1)+1);
RC0=zeros(ynum(1)-1,xnum(1)-1);
% Grid points cycle
for i=1:1:ynum(1);
    for j=1:1:xnum(1);
        % Relative distances for (i,j) basic node inside the grid
        dx=(j-1)*xstp(1)/xsize;
        dy=(i-1)*ystp(1)/ysize;
        % Density, viscosity structure (falling block)
        rho1(i,j)=rhomedium;
        etas1(i,j)=etamedium;
        if(dx>=0.4 && dx<=0.6 && dy>=0.4 && dy<=0.6)
            rho1(i,j)=rhoblock;
            etas1(i,j)=etablock;
        end
        % Viscosity for normal stress (etan) for cells (in pressure nodes)
        if (i<ynum(1) && j<xnum(1))
            % Relative distances for (i,j) pressure nodes inside the grid
            dx=(j-0.5)*xstp(1)/xsize;
            dy=(i-0.5)*ystp(1)/ysize;
            etan1(i,j)=etamedium;
            if(dx>=0.4 && dx<=0.6 && dy>=0.4 && dy<=0.6)
                etan1(i,j)=etablock;
            end
        end
        % Pressure update in internal nodes (lithostatic pressure)
        if(i<ynum(1) && j<xnum(1))
            if(i==1)
                pr0(i,j)=prfirst;
            else
                pr0(i,j)=pr0(i-1,j)+ystp(1)*g*(rho1(i,j)+rho1(i,j+1))/2;
            end
        end
        % Right part of x-Stokes Equation
        if(j>1 && i>1 && j<xnum(1))
            RX0(i,j)=0;
        end
        % Right part of y-Stokes Equation
        if(j>1 && i>1 && i<ynum(1))
            RY0(i,j)=-g*(rho1(i,j)+rho1(i,j-1))/2;
        end
    end
end

%Computing, Saving initial state of residuals
% by calling function Stokes_smoother()
[vx0,resx0,vy0,resy0,pr0,resc0]=Stokes_Continuity_viscous_smoother(prfirst,etas1,etan1,0,krelaxs,krelaxc,xnum(k),ynum(k),xstp(k),ystp(k),RX0,vx0,RY0,vy0,RC0,pr0);
RX1=resx0;
RY1=resy0;
RC1=resc0;

 
% Interpolation (restriction) of viscosity to coarser levels
for k=1:1:levelnum-1;
    % Action depends on the level of resolution
    switch k
            
        % Level 1 (principal grid)
        case 1
        if (levelnum>k)
            [etas2 etan2]=Viscosity_restriction(etamean,k,xnum,ynum,xstp,ystp,etas1,etan1);
        end
        % Level 2 
        case 2
        if (levelnum>k)
            [etas3 etan3]=Viscosity_restriction(etamean,k,xnum,ynum,xstp,ystp,etas2,etan2);
        end
        % Level 3 
        case 3
        if (levelnum>k)
            [etas4 etan4]=Viscosity_restriction(etamean,k,xnum,ynum,xstp,ystp,etas3,etan3);
        end
        % Level 4 
        case 4
        if (levelnum>k)
            [etas5 etan5]=Viscosity_restriction(etamean,k,xnum,ynum,xstp,ystp,etas4,etan4);
        end
        % Level 5 
        case 5
        if (levelnum>k)
            [etas6 etan6]=Viscosity_restriction(etamean,k,xnum,ynum,xstp,ystp,etas5,etan5);
        end
        % Level 6 
        case 6
        if (levelnum>k)
            [etas7 etan7]=Viscosity_restriction(etamean,k,xnum,ynum,xstp,ystp,etas6,etan6);
        end
        % Level 7 
        case 7
        if (levelnum>k)
            [etas8 etan8]=Viscosity_restriction(etamean,k,xnum,ynum,xstp,ystp,etas7,etan7);
        end
        % No further restriction from last Level 8
    end
end


% Saving viscosity
if (levelnum>0) 
    etaso1=etas1; 
    etano1=etan1; 
end
if (levelnum>1) 
    etaso2=etas2; 
    etano2=etan2; 
end
if (levelnum>2) 
    etaso3=etas3; 
    etano3=etan3; 
end
if (levelnum>3) 
    etaso4=etas4; 
    etano4=etan4; 
end
if (levelnum>4) 
    etaso5=etas5; 
    etano5=etan5; 
end
if (levelnum>5) 
    etaso6=etas6; 
    etano6=etan6; 
end
if (levelnum>6) 
    etaso7=etas7; 
    etano7=etan7; 
end
if (levelnum>7) 
    etaso8=etas8; 
    etano8=etan8; 
end


% Renormalizing viscosity
if (etamaxcur<etamax || etamincur>etamin)
    % Normalization coefficient
    etakf=log(etamaxcur/etamincur)/log(etamax/etamin);
    % Recomputing viscosity
    if (levelnum>0) 
        etas1=etamincur*exp(log(etaso1/etamin)*etakf);
        etan1=etamincur*exp(log(etano1/etamin)*etakf);
    end
    if (levelnum>1) 
        etas2=etamincur*exp(log(etaso2/etamin)*etakf);
        etan2=etamincur*exp(log(etano2/etamin)*etakf);
    end
    if (levelnum>2) 
        etas3=etamincur*exp(log(etaso3/etamin)*etakf);
        etan3=etamincur*exp(log(etano3/etamin)*etakf);
    end
    if (levelnum>3) 
        etas4=etamincur*exp(log(etaso4/etamin)*etakf);
        etan4=etamincur*exp(log(etano4/etamin)*etakf);
    end
    if (levelnum>4) 
        etas5=etamincur*exp(log(etaso5/etamin)*etakf);
        etan5=etamincur*exp(log(etano5/etamin)*etakf);
    end
    if (levelnum>5) 
        etas6=etamincur*exp(log(etaso6/etamin)*etakf);
        etan6=etamincur*exp(log(etano6/etamin)*etakf);
    end
    if (levelnum>6) 
        etas7=etamincur*exp(log(etaso7/etamin)*etakf);
        etan7=etamincur*exp(log(etano7/etamin)*etakf);
    end
    if (levelnum>7) 
        etas8=etamincur*exp(log(etaso8/etamin)*etakf);
        etan8=etamincur*exp(log(etano8/etamin)*etakf);
    end
end


% Main Multigrid cycle
% We put 8 levels (for the case if we want to increase resolution)
% We use different arrays for different levels 
% since grid resolutions is different
% Multigrid schedule = V-cycle
niterbig=1;
ynreset=0
for niter=1:1:inum;
    
    
    % Smoothing+restriction cycle
    for k=1:1:levelnum;
        % Action depends on the level of resolution
        switch k
            
            % Leve 1 (principal grid)
            case 1
            % Smoothing operation: solving of Stokes and Continuity equations on nodes
            % and computing residuals
            % by calling function Stokes_smoother()
            [vx1,resx1,vy1,resy1,pr1,resc1]=Stokes_Continuity_viscous_smoother(0,etas1,etan1,iternum(k),krelaxs,krelaxc,xnum(k),ynum(k),xstp(k),ystp(k),RX1,vx1,RY1,vy1,RC1,pr1);
            % Restriction operation: 
            % Interpolating residuals to coarcer level (k+1) 
            % to produce right parts for this coarcer level
            if (levelnum>k)
                [RX2,RY2,RC2]=Stokes_Continuity_viscous_restriction(k,xnum,ynum,xstp,ystp,resx1,resy1,resc1,etan1,etan2);
            end
        
            % Level 2
            case 2
            %Making initial approximation for vx,vy, P corrections (zeros)
            vx2=zeros(ynum(k)+1,xnum(k));
            vy2=zeros(ynum(k),xnum(k)+1);
            pr2=zeros(ynum(k)-1,xnum(k)-1);
            % Smoothing operation: 
            [vx2,resx2,vy2,resy2,pr2,resc2]=Stokes_Continuity_viscous_smoother(0,etas2,etan2,iternum(k),krelaxs,krelaxc,xnum(k),ynum(k),xstp(k),ystp(k),RX2,vx2,RY2,vy2,RC2,pr2);
            % Restriction operation: 
            if (levelnum>k)
                [RX3,RY3,RC3]=Stokes_Continuity_viscous_restriction(k,xnum,ynum,xstp,ystp,resx2,resy2,resc2,etan2,etan3);
            end
        
            % Level 3
            case 3
            %Making initial approximation for vx,vy, P corrections (zeros)
            vx3=zeros(ynum(k)+1,xnum(k));
            vy3=zeros(ynum(k),xnum(k)+1);
            pr3=zeros(ynum(k)-1,xnum(k)-1);
            % Smoothing operation: 
            [vx3,resx3,vy3,resy3,pr3,resc3]=Stokes_Continuity_viscous_smoother(0,etas3,etan3,iternum(k),krelaxs,krelaxc,xnum(k),ynum(k),xstp(k),ystp(k),RX3,vx3,RY3,vy3,RC3,pr3);
            % Restriction operation: 
            if (levelnum>k)
                [RX4,RY4,RC4]=Stokes_Continuity_viscous_restriction(k,xnum,ynum,xstp,ystp,resx3,resy3,resc3,etan3,etan4);
            end
            
            % Level 4
            case 4
            %Making initial approximation for vx,vy, P corrections (zeros)
            vx4=zeros(ynum(k)+1,xnum(k));
            vy4=zeros(ynum(k),xnum(k)+1);
            pr4=zeros(ynum(k)-1,xnum(k)-1);
            % Smoothing operation: 
            [vx4,resx4,vy4,resy4,pr4,resc4]=Stokes_Continuity_viscous_smoother(0,etas4,etan4,iternum(k),krelaxs,krelaxc,xnum(k),ynum(k),xstp(k),ystp(k),RX4,vx4,RY4,vy4,RC4,pr4);
            % Restriction operation: 
            if (levelnum>k)
                [RX5,RY5,RC5]=Stokes_Continuity_viscous_restriction(k,xnum,ynum,xstp,ystp,resx4,resy4,resc4,etan4,etan5);
            end
            
             % Level 5
            case 5
            %Making initial approximation for vx,vy, P corrections (zeros)
            vx5=zeros(ynum(k)+1,xnum(k));
            vy5=zeros(ynum(k),xnum(k)+1);
            pr5=zeros(ynum(k)-1,xnum(k)-1);
            % Smoothing operation: 
            [vx5,resx5,vy5,resy5,pr5,resc5]=Stokes_Continuity_viscous_smoother(0,etas5,etan5,iternum(k),krelaxs,krelaxc,xnum(k),ynum(k),xstp(k),ystp(k),RX5,vx5,RY5,vy5,RC5,pr5);
            % Restriction operation: 
            if (levelnum>k)
                [RX6,RY6,RC6]=Stokes_Continuity_viscous_restriction(k,xnum,ynum,xstp,ystp,resx5,resy5,resc5,etan5,etan6);
            end
             
            % Level 6
            case 6
            %Making initial approximation for vx,vy, P corrections (zeros)
            vx6=zeros(ynum(k)+1,xnum(k));
            vy6=zeros(ynum(k),xnum(k)+1);
            pr6=zeros(ynum(k)-1,xnum(k)-1);
            % Smoothing operation: 
            [vx6,resx6,vy6,resy6,pr6,resc6]=Stokes_Continuity_viscous_smoother(0,etas6,etan6,iternum(k),krelaxs,krelaxc,xnum(k),ynum(k),xstp(k),ystp(k),RX6,vx6,RY6,vy6,RC6,pr6);
            % Restriction operation: 
            if (levelnum>k)
                [RX7,RY7,RC7]=Stokes_Continuity_viscous_restriction(k,xnum,ynum,xstp,ystp,resx6,resy6,resc6,etan6,etan7);
            end
             
            % Level 7
            case 7
            %Making initial approximation for vx,vy, P corrections (zeros)
            vx7=zeros(ynum(k)+1,xnum(k));
            vy7=zeros(ynum(k),xnum(k)+1);
            pr7=zeros(ynum(k)-1,xnum(k)-1);
            % Smoothing operation: 
            [vx7,resx7,vy7,resy7,pr7,resc7]=Stokes_Continuity_viscous_smoother(0,etas7,etan7,iternum(k),krelaxs,krelaxc,xnum(k),ynum(k),xstp(k),ystp(k),RX7,vx7,RY7,vy7,RC7,pr7);
            % Restriction operation: 
            if (levelnum>k)
                [RX8,RY8,RC8]=Stokes_Continuity_viscous_restriction(k,xnum,ynum,xstp,ystp,resx7,resy7,resc7,etan7,etan8);
            end
             
            % Level 8
            case 8
            %Making initial approximation for vx,vy, P corrections (zeros)
            vx8=zeros(ynum(k)+1,xnum(k));
            vy8=zeros(ynum(k),xnum(k)+1);
            pr8=zeros(ynum(k)-1,xnum(k)-1);
            % Smoothing operation: 
            [vx8,resx8,vy8,resy8,pr8,resc8]=Stokes_Continuity_viscous_smoother(0,etas8,etan8,iternum(k),krelaxs,krelaxc,xnum(k),ynum(k),xstp(k),ystp(k),RX8,vx8,RY8,vy8,RC8,pr8);
            % No Restriction operation on this last level 

        end
    end
    
    % Smoothing+prolongation cycle
    for k=levelnum:-1:1;
        % Action depends on the level of resolution
        switch k

            % Leve 1 (principal grid)
            case 1
            % Smoothing operation: solving of Stokes and Continuity equations on nodes
            % and computing residuals
            % by calling function Stokes_smoother()
            [vx1,resx1,vy1,resy1,pr1,resc1]=Stokes_Continuity_viscous_smoother(0,etas1,etan1,iternum(k),krelaxs,krelaxc,xnum(k),ynum(k),xstp(k),ystp(k),RX1,vx1,RY1,vy1,RC1,pr1);
            % No prolongation operation on finest level
        
            % Level 2
            case 2
            % Smoothing operation: 
            [vx2,resx2,vy2,resy2,pr2,resc2]=Stokes_Continuity_viscous_smoother(0,etas2,etan2,iternum(k),krelaxs,krelaxc,xnum(k),ynum(k),xstp(k),ystp(k),RX2,vx2,RY2,vy2,RC2,pr2);
            % Prolongation operation: 
            % Interpolating dvx,dvy,dpr corrections to finer level (k+1) 
            % and computing solution update for this finer level
            [dvx1,dvy1,dpr1]=Stokes_Continuity_prolongation(k,xnum,ynum,xstp,ystp,vx2,vy2,pr2);
            vx1=vx1+dvx1*vkoef;
            vy1=vy1+dvy1*vkoef;
            pr1=pr1+dpr1*pkoef;
        
            % Level 3
            case 3
            % Smoothing operation: 
            [vx3,resx3,vy3,resy3,pr3,resc3]=Stokes_Continuity_viscous_smoother(0,etas3,etan3,iternum(k),krelaxs,krelaxc,xnum(k),ynum(k),xstp(k),ystp(k),RX3,vx3,RY3,vy3,RC3,pr3);
            % Prolongation operation: 
            [dvx2,dvy2,dpr2]=Stokes_Continuity_prolongation(k,xnum,ynum,xstp,ystp,vx3,vy3,pr3);
            vx2=vx2+dvx2*vkoef;
            vy2=vy2+dvy2*vkoef;
            pr2=pr2+dpr2*pkoef;
            
            % Level 4
            case 4
            % Smoothing operation: 
            [vx4,resx4,vy4,resy4,pr4,resc4]=Stokes_Continuity_viscous_smoother(0,etas4,etan4,iternum(k),krelaxs,krelaxc,xnum(k),ynum(k),xstp(k),ystp(k),RX4,vx4,RY4,vy4,RC4,pr4);
            % Prolongation operation: 
            [dvx3,dvy3,dpr3]=Stokes_Continuity_prolongation(k,xnum,ynum,xstp,ystp,vx4,vy4,pr4);
            vx3=vx3+dvx3*vkoef;
            vy3=vy3+dvy3*vkoef;
            pr3=pr3+dpr3*pkoef;
            
             % Level 5
            case 5
            % Smoothing operation: 
            [vx5,resx5,vy5,resy5,pr5,resc5]=Stokes_Continuity_viscous_smoother(0,etas5,etan5,iternum(k),krelaxs,krelaxc,xnum(k),ynum(k),xstp(k),ystp(k),RX5,vx5,RY5,vy5,RC5,pr5);
            % Prolongation operation: 
            [dvx4,dvy4,dpr4]=Stokes_Continuity_prolongation(k,xnum,ynum,xstp,ystp,vx5,vy5,pr5);
            vx4=vx4+dvx4*vkoef;
            vy4=vy4+dvy4*vkoef;
            pr4=pr4+dpr4*pkoef;
            
            
            % Level 6
            case 6
            % Smoothing operation: 
            [vx6,resx6,vy6,resy6,pr6,resc6]=Stokes_Continuity_viscous_smoother(0,etas6,etan6,iternum(k),krelaxs,krelaxc,xnum(k),ynum(k),xstp(k),ystp(k),RX6,vx6,RY6,vy6,RC6,pr6);
            % Prolongation operation: 
            [dvx5,dvy5,dpr5]=Stokes_Continuity_prolongation(k,xnum,ynum,xstp,ystp,vx6,vy6,pr6);
            vx5=vx5+dvx5*vkoef;
            vy5=vy5+dvy5*vkoef;
            pr5=pr5+dpr5*pkoef;
            
            % Level 7
            case 7
            % Smoothing operation: 
            [vx7,resx7,vy7,resy7,pr7,resc7]=Stokes_Continuity_viscous_smoother(0,etas7,etan7,iternum(k),krelaxs,krelaxc,xnum(k),ynum(k),xstp(k),ystp(k),RX7,vx7,RY7,vy7,RC7,pr7);
            % Prolongation operation: 
            [dvx6,dvy6,dpr6]=Stokes_Continuity_prolongation(k,xnum,ynum,xstp,ystp,vx7,vy7,pr7);
            vx6=vx6+dvx6*vkoef;
            vy6=vy6+dvy6*vkoef;
            pr6=pr6+dpr6*pkoef;
            
            % Level 8
            case 8
            % Smoothing operation: 
            [vx8,resx8,vy8,resy8,pr8,resc8]=Stokes_Continuity_viscous_smoother(0,etas8,etan8,iternum(k),krelaxs,krelaxc,xnum(k),ynum(k),xstp(k),ystp(k),RX8,vx8,RY8,vy8,RC8,pr8);
            % Prolongation operation: 
            [dvx7,dvy7,dpr7]=Stokes_Continuity_prolongation(k,xnum,ynum,xstp,ystp,vx8,vy8,pr8);
            vx7=vx7+dvx7*vkoef;
            vy7=vy7+dvy7*vkoef;
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
    resx01=resx1/stokesscale;
    surf(resx01);
    shading interp;
    light;
    lighting phong;
    axis tight;
    zlabel('residual x-Stokes')
%     title(['x-Stokes residual, Gauss-Seidel cycle = ',num2str(niter)])
%     colorbar;
    
    % Plotting Residuals for y-Stokes as surface
    subplot(2,3,2)
    resy01=resy1/stokesscale;
    surf(resy01);
    shading interp;
    light;
    lighting phong;
    axis tight;
    zlabel('residual Y-stokes')
%     title(['y-Stokes residual, Gauss-Seidel cycle = ',num2str(niter)])
%     colorbar;
    
    % Plotting Residuals for Continuity as surface
    subplot(2,3,3)
    resc01=resc1/continscale;
    surf(resc01);
    shading interp;
    light;
    lighting phong;
    axis tight;
    zlabel('residual continuity')
%     title(['Continuity residual, Gauss-Seidel cycle = ',num2str(niter)])
%     colorbar;
    
    % Plotting Pressure
    subplot(2,3,4)
    surf(vx1);
    shading interp;
    light;
    lighting phong;
    axis tight;
    zlabel(['Vx, cycle ',num2str(niter)])
%     title(['P, Multigrid V-cycle = ',num2str(niter)])
%     colorbar;

    % Plotting vx
    subplot(2,3,5)
    surf(vy1);
    shading interp;
    light;
    lighting phong;
    axis tight;
    zlabel(['Vy, cycle ',num2str(niter)])
%     title(['Pressure, Multigrid V-cycle = ',num2str(niter)])
%     colorbar;

    % Computing mean square residuals
    resx001(niter)=0;
    resy001(niter)=0;
    resc001(niter)=0;
    for i=1:1:ynum(1)
        for j=1:1:xnum(1)
            % x-Stokes
            if (i>1 && j>1 && j<xnum(1))
                resx001(niter)=resx001(niter)+resx01(i,j)^2;
            end
            % y-Stokes
            if (i>1 && j>1 && i<ynum(1))
                resy001(niter)=resy001(niter)+resy01(i,j)^2;
            end
            % Continuity
            if (j<xnum(1) && i<ynum(1))
                resc001(niter)=resc001(niter)+resc01(i,j)^2;
            end
        end
    end
    resx001(niter)=log10((resx001(niter)/((ynum(1)-1)*(xnum(1)-2)))^0.5);
    resy001(niter)=log10((resy001(niter)/((ynum(1)-2)*(xnum(1)-1)))^0.5);
    resc001(niter)=log10((resc001(niter)/((ynum(1)-1)*(xnum(1)-1)))^0.5);
   

    % Plotting Mean residuals
    subplot(2,3,6)
    plot(resc001, 'k');
    hold on
    plot(resx001, 'b');
    hold on
    plot(resy001, 'r');
    xlabel('V-cycles')
    ylabel('log(Residuals)')

    zlabel(['RMS errors, cycle ',num2str(niter)])

    % Delay program execution
    pause(0.1)
    
    % Renormalizing viscosity
    % Add cycle counter
    etannn=etannn+1;
    if(etannn>etacyc)
        % Reset counter
        etannn=1;
        % Change lower viscosity limit
        etamaxcur=etamaxcur*etamaxstp;
        if (etamaxcur>=etamax)
            ynreset=ynreset+1;
        end
        if (etamaxcur>etamax)
            etamaxcur=etamax;
        end
        if (ynreset>1)
            %Viscosity reset
            etamaxcur=etamaxbeg;
            ynreset=0;
            %Computing, Saving initial state of residuals
            % by calling function Stokes_smoother()
            vx0=vx0+vx1;
            vy0=vy0+vy1;
            pr0=pr0+pr1;
            [vx0,resx0,vy0,resy0,pr0,resc0]=Stokes_Continuity_viscous_smoother(prfirst,etaso1,etano1,0,krelaxs,krelaxc,xnum(k),ynum(k),xstp(k),ystp(k),RX0,vx0,RY0,vy0,RC0,pr0);
            RX1=resx0;
            RY1=resy0;
            RC1=resc0;   
            %Making initial approximation for vx,vy, P corrections (zeros)
            vx1=zeros(ynum(1)+1,xnum(1));
            vy1=zeros(ynum(1),xnum(1)+1);
            pr1=zeros(ynum(1)-1,xnum(1)-1);
    % Defining scale for Stokes residuals from y-Stokes equation
    % dSIGMAij/dj-dP/di=-RHO*gi=0  => Stokes scale=abs(RHO*gi)
    stokesscale= rhoblock*g;
    % Defining scale for Continuity residuals from y-Stokes equation
    % dvx/dx+dvy/dy=0 can be transformed to 2ETA(dvx/dx+dvy/dy)/dx=0 
    % which is similar to dSIGMAij/dj and has scale given above
    % therefore continuity scale = scale=abs(RHO*gi/ETA*dx)
    continscale= rhoblock*g/etamedium*xstp(1);


    % Plotting Residuals for x-Stokes as surface
    figure(2);
    subplot(2,3,1)
    resx01=resx0/stokesscale;
    surf(resx01);
    shading interp;
    light;
    lighting phong;
    axis tight;
    zlabel('residual x-Stokes')
%     title(['x-Stokes residual, Gauss-Seidel cycle = ',num2str(niter)])
%     colorbar;
    
    % Plotting Residuals for y-Stokes as surface
    subplot(2,3,2)
    resy01=resy0/stokesscale;
    surf(resy01);
    shading interp;
    light;
    lighting phong;
    axis tight;
    zlabel('residual Y-stokes')
%     title(['y-Stokes residual, Gauss-Seidel cycle = ',num2str(niter)])
%     colorbar;
    
    % Plotting Residuals for Continuity as surface
    subplot(2,3,3)
    resc01=resc0/continscale;
    surf(resc01);
    shading interp;
    light;
    lighting phong;
    axis tight;
    zlabel('residual continuity')
%     title(['Continuity residual, Gauss-Seidel cycle = ',num2str(niter)])
%     colorbar;
    
    % Plotting Pressure
    subplot(2,3,4)
    surf(vx0);
    shading interp;
    light;
    lighting phong;
    axis tight;
    zlabel(['Vx, cycle ',num2str(niter)])
%     title(['P, Multigrid V-cycle = ',num2str(niter)])
%     colorbar;

    % Plotting vx
    subplot(2,3,5)
    surf(vy0);
    shading interp;
    light;
    lighting phong;
    axis tight;
    zlabel(['Vy, cycle ',num2str(niter)])
%     title(['Pressure, Multigrid V-cycle = ',num2str(niter)])
%     colorbar;

    % Computing mean square residuals
    resx000(niterbig)=0;
    resy000(niterbig)=0;
    resc000(niterbig)=0;
    for i=1:1:ynum(1)
        for j=1:1:xnum(1)
            % x-Stokes
            if (i>1 && j>1 && j<xnum(1))
                resx000(niterbig)=resx000(niterbig)+resx01(i,j)^2;
            end
            % y-Stokes
            if (i>1 && j>1 && i<ynum(1))
                resy000(niterbig)=resy000(niterbig)+resy01(i,j)^2;
            end
            % Continuity
            if (j<xnum(1) && i<ynum(1))
                resc000(niterbig)=resc000(niterbig)+resc01(i,j)^2;
            end
        end
    end
    resx000(niterbig)=log10((resx000(niterbig)/((ynum(1)-1)*(xnum(1)-2)))^0.5);
    resy000(niterbig)=log10((resy000(niterbig)/((ynum(1)-2)*(xnum(1)-1)))^0.5);
    resc000(niterbig)=log10((resc000(niterbig)/((ynum(1)-1)*(xnum(1)-1)))^0.5);
        % Plotting Mean residuals
    subplot(2,3,6)
    plot(resc000, 'k');
    hold on
    plot(resx000, 'b');
    hold on
    plot(resy000, 'r');
    xlabel('V-cycles')
    ylabel('log(Residuals)')

    %Update counter
    niterbig=niterbig+1;
    pause(3);
        end
        
        % Normalization coefficient
        etakf=log(etamaxcur/etamincur)/log(etamax/etamin);
        % Recomputing viscosity
        if (levelnum>0) 
            etas1=etamincur*exp(log(etaso1/etamin)*etakf);
            etan1=etamincur*exp(log(etano1/etamin)*etakf);
        end
        if (levelnum>1) 
            etas2=etamincur*exp(log(etaso2/etamin)*etakf);
            etan2=etamincur*exp(log(etano2/etamin)*etakf);
        end
        if (levelnum>2) 
            etas3=etamincur*exp(log(etaso3/etamin)*etakf);
            etan3=etamincur*exp(log(etano3/etamin)*etakf);
        end
        if (levelnum>3) 
            etas4=etamincur*exp(log(etaso4/etamin)*etakf);
            etan4=etamincur*exp(log(etano4/etamin)*etakf);
        end
        if (levelnum>4) 
            etas5=etamincur*exp(log(etaso5/etamin)*etakf);
            etan5=etamincur*exp(log(etano5/etamin)*etakf);
        end
        if (levelnum>5) 
            etas6=etamincur*exp(log(etaso6/etamin)*etakf);
            etan6=etamincur*exp(log(etano6/etamin)*etakf);
        end
        if (levelnum>6) 
            etas7=etamincur*exp(log(etaso7/etamin)*etakf);
            etan7=etamincur*exp(log(etano7/etamin)*etakf);
        end
        if (levelnum>7) 
            etas8=etamincur*exp(log(etaso8/etamin)*etakf);
            etan8=etamincur*exp(log(etano8/etamin)*etakf);
        end
    end


end

