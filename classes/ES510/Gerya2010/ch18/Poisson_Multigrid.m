% Solution of 2D Poisson equation 
% with Multigrid based on V-cycle
% using external functions Poisson_smoother.m, Poisson_restriction.m,
% Poisson_prolongation.m
% resolution between multigrid levels changes by factor of two

% Clearing all variables and arrays
clear all;
% Clearing figures
clf;

% Model parameters
% Radius of the planet
r=6000000.0;
% Density of the planet
rhoplanet=6000.0;
% Gravity constant
G=6.67e-11;
% Model size, m
xsize=18000000.0;
ysize=18000000.0;

% Multigrid parameters
% Number of resolution levels
levelnum=6;
% Resolution (number of nodal points) on the coarsest (last) grid
xnumcur=6+1;
ynumcur=6+1;

% Iteration Parameters
% Total number of iteration cycles
inum=12;
% Relaxation coefficient for Gauss-Seidel iterations
krelax=1.75;
% Number of smoothing iterations
iternum(1)=5;
iternum(2)=5*2;
iternum(3)=5*4;
iternum(4)=5*8;
iternum(5)=5*16;
iternum(6)=5*32;
iternum(7)=5*64;
iternum(8)=5*128;


% Defining number of iterations, resolutions and gridsteps on all levels
% Number 
for k=levelnum:-1:1
    % Numerical resolution for current level
    xnum(k)=xnumcur;
    ynum(k)=ynumcur;
    % Grid steps
    xstp(k)=xsize./(xnum(k)-1);
    ystp(k)=ysize./(ynum(k)-1);
    % Numerical resolution for the next finer level
    % Doubling number of CELLS
    xnumcur=(xnumcur-1)*2+1;
    ynumcur=(ynumcur-1)*2+1;
end

% FINEST (PRINCIPAL) GRID
% Defining density structure rho()
% Defining initial guesses for gravity potential phi()
% Computing right part of the Poisson equation R()
% Grid points cycle
for i=1:1:ynum(1);
    for j=1:1:xnum(1);
        % Gravity potential
        phi1(i,j)=0;
        % Check distance of (i,j) node from the grid center
        dx=(j-1).*xstp(1)-xsize./2;
        dy=(i-1).*ystp(1)-ysize./2;
        dist=(dx.*dx+dy.*dy)^0.5;
        % Density 
        rho(i,j)=0;
        if (dist<r)
            rho(i,j)=rhoplanet;
        end
        % Right part of Poisson equation
        R1(i,j)= 4.0.*pi.*G.*rho(i,j);
    end
end

% Figures counter
fignum=1;

% Main Multigrid cycle
% We put 8 levels (for the case if we want to increase resolution)
% We use different arrays for different levels 
% since grid resolutions is different
% Multigrid schedule = V-cycle 
for niter=1:1:inum;
    
    % Smoothing+restriction cycle
    for k=1:1:levelnum;
        % Action depends on the level of resolution
        switch k
            
            % Leve 1 (principal grid)
            case 1
            % Smoothing operation: solving of Poisson equation on nodes
            % and computing residuals
            % by calling function Poisson_smoother()
            [phi1,residual1]=Poisson_smoother(iternum(k),krelax,xnum(k),ynum(k),xstp(k),ystp(k),R1,phi1);
            % Restriction operation: 
            % Interpolating residuals to coarcer level (k+1) 
            % to produce right parts for this coarcer level
            if (levelnum>k)
                [R2]=Poisson_restriction(k,xnum,ynum,xstp,ystp,residual1);
            end
        
            % Level 2
            case 2
            %Making initial approximation for phi corrections (zeros)
            phi2=zeros(ynum(k),xnum(k));
            % Smoothing operation: 
            [phi2,residual2]=Poisson_smoother(iternum(k),krelax,xnum(k),ynum(k),xstp(k),ystp(k),R2,phi2);
            % Restriction operation: 
            if (levelnum>k)
                [R3]=Poisson_restriction(k,xnum,ynum,xstp,ystp,residual2);
            end
        
            % Level 3
            case 3
            %Making initial approximation for phi corrections (zeros)
            phi3=zeros(ynum(k),xnum(k));
            % Smoothing operation: 
            [phi3,residual3]=Poisson_smoother(iternum(k),krelax,xnum(k),ynum(k),xstp(k),ystp(k),R3,phi3);
            % Restriction operation: 
            if (levelnum>k)
                [R4]=Poisson_restriction(k,xnum,ynum,xstp,ystp,residual3);
            end
            
            % Level 4
            case 4
            %Making initial approximation for phi corrections (zeros)
            phi4=zeros(ynum(k),xnum(k));
            % Smoothing operation: 
            [phi4,residual4]=Poisson_smoother(iternum(k),krelax,xnum(k),ynum(k),xstp(k),ystp(k),R4,phi4);
            % Restriction operation: 
            if (levelnum>k)
                [R5]=Poisson_restriction(k,xnum,ynum,xstp,ystp,residual4);
            end
            
             % Level 5
            case 5
            %Making initial approximation for phi corrections (zeros)
            phi5=zeros(ynum(k),xnum(k));
            % Smoothing operation: 
            [phi5,residual5]=Poisson_smoother(iternum(k),krelax,xnum(k),ynum(k),xstp(k),ystp(k),R5,phi5);
            % Restriction operation: 
            if (levelnum>k)
                [R6]=Poisson_restriction(k,xnum,ynum,xstp,ystp,residual5);
            end
             
            % Level 6
            case 6
            %Making initial approximation for phi corrections (zeros)
            phi6=zeros(ynum(k),xnum(k));
            % Smoothing operation: 
            [phi6,residual6]=Poisson_smoother(iternum(k),krelax,xnum(k),ynum(k),xstp(k),ystp(k),R6,phi6);
              % Restriction operation: 
            if (levelnum>k)
                [R7]=Poisson_restriction(k,xnum,ynum,xstp,ystp,residual6);
            end
             
            % Level 7
            case 7
            %Making initial approximation for phi corrections (zeros)
            phi7=zeros(ynum(k),xnum(k));
            % Smoothing operation: 
            [phi7,residual7]=Poisson_smoother(iternum(k),krelax,xnum(k),ynum(k),xstp(k),ystp(k),R7,phi7);
              % Restriction operation: 
            if (levelnum>k)
                [R8]=Poisson_restriction(k,xnum,ynum,xstp,ystp,residual7);
            end
             
            % Level 8
            case 8
            %Making initial approximation for phi corrections (zeros)
            phi8=zeros(ynum(k),xnum(k));
            % Smoothing operation: 
            [phi8,residual8]=Poisson_smoother(iternum(k),krelax,xnum(k),ynum(k),xstp(k),ystp(k),R8,phi8);
            % No Restriction operation on this last level 

        end
    end
   
    % Smoothing+prolongation cycle
    for k=levelnum:-1:1;
        % Action depends on the level of resolution
        switch k

            % Leve 1 (principal grid)
            case 1
            % Smoothing operation: solving of Poisson equation on nodes
            % and computing residuals
            % by calling function Poisson_smoother()
            [phi1,residual1]=Poisson_smoother(iternum(k),krelax,xnum(k),ynum(k),xstp(k),ystp(k),R1,phi1);
            % No prolongation operation on finest level
        
            % Level 2
            case 2
            % Smoothing operation: 
            [phi2,residual2]=Poisson_smoother(iternum(k),krelax,xnum(k),ynum(k),xstp(k),ystp(k),R2,phi2);
            % Prolongation operation: 
            % Interpolating phi corrections to finer level (k+1) 
            % and computing solution update for this finer level
            [dphi1]=Poisson_prolongation(k,xnum,ynum,xstp,ystp,phi2);
            phi1=phi1+dphi1;
        
            % Level 3
            case 3
            % Smoothing operation: 
            [phi3,residual3]=Poisson_smoother(iternum(k),krelax,xnum(k),ynum(k),xstp(k),ystp(k),R3,phi3);
            % Prolongation operation: 
            [dphi2]=Poisson_prolongation(k,xnum,ynum,xstp,ystp,phi3);
            phi2=phi2+dphi2;
            
            % Level 4
            case 4
            % Smoothing operation: 
            [phi4,residual4]=Poisson_smoother(iternum(k),krelax,xnum(k),ynum(k),xstp(k),ystp(k),R4,phi4);
            % Prolongation operation: 
            [dphi3]=Poisson_prolongation(k,xnum,ynum,xstp,ystp,phi4);
            phi3=phi3+dphi3;
            
             % Level 5
            case 5
            % Smoothing operation: 
            [phi5,residual5]=Poisson_smoother(iternum(k),krelax,xnum(k),ynum(k),xstp(k),ystp(k),R5,phi5);
            % Prolongation operation: 
            [dphi4]=Poisson_prolongation(k,xnum,ynum,xstp,ystp,phi5);
            phi4=phi4+dphi4;
            
            % Level 6
            case 6
            % Smoothing operation: 
            [phi6,residual6]=Poisson_smoother(iternum(k),krelax,xnum(k),ynum(k),xstp(k),ystp(k),R6,phi6);
            % Prolongation operation: 
            [dphi5]=Poisson_prolongation(k,xnum,ynum,xstp,ystp,phi6);
            phi5=phi5+dphi5;
            
            % Level 7
            case 7
            % Smoothing operation: 
            [phi7,residual7]=Poisson_smoother(iternum(k),krelax,xnum(k),ynum(k),xstp(k),ystp(k),R7,phi7);
            % Prolongation operation: 
            [dphi6]=Poisson_prolongation(k,xnum,ynum,xstp,ystp,phi7);
            phi6=phi6+dphi6;
            
            % Level 8
            case 8
            % Smoothing operation: 
            [phi8,residual8]=Poisson_smoother(iternum(k),krelax,xnum(k),ynum(k),xstp(k),ystp(k),R8,phi8);
            % Prolongation operation: 
            [dphi7]=Poisson_prolongation(k,xnum,ynum,xstp,ystp,phi8);
            phi7=phi7+dphi7;
        end
    end 
 
    % Plotting Residuals for selected stages
    n=niter;
    if(n==1 || n==3 || n==6 || n==12)
        % Normalizing residual
        kfr=4.0.*pi.*G*rhoplanet;
        residual0=residual1/kfr;
        % Plotting normalized residuals as surface
        figure(fignum);
        fignum=fignum+1;
        surf(residual0);
        colormap(Gray);
        shading interp;
        light;
        lighting phong;
        axis tight;
        zlabel('residual');
        title(['Solution of Poisson equation, iteration = ',num2str(n*iternum(1))]);
        pause(1);
    end

end

