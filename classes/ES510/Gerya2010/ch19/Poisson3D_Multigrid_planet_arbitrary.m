% Solution of Poisson equation in 3D with Multigrid
% based on V-cycle
% by using external function Poisson3D_smoother()
% Density distribution in the model corresponds to a planet
% grids resolution is arbitryry chosen

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
zsize=18000000.0;
% Radius for gravity potential boundary
gradius=0.9999*xsize/2;

% Multigrid parameters
% Number of resolution levels
levelnum=4;

% Iteration Parameters
% Total number of iteration cycles
inum=30;
% Relaxation coefficient for Gauss-Seidel iterations
krelax=1.5;
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

xnum(1)=49;
ynum(1)=49;
znum(1)=49;
xnum(2)=25;
ynum(2)=25;
znum(2)=25;
xnum(3)=13;
ynum(3)=13;
znum(3)=13;
xnum(4)=7;
ynum(4)=7;
znum(4)=7;


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
% Defining initial guesses for gravity potential phi()
% Computing right part of the Poisson equation R()
% Grid points cycle
for i=1:1:ynum(1);
    for j=1:1:xnum(1);
        for k=1:1:znum(1);

            % Gravity potential
            phi1(i,j,k)=0;
            % Check distance of (i,j) node from the grid center
            dx=(j-1).*xstp(1)-xsize./2;
            dy=(i-1).*ystp(1)-ysize./2;
            dz=(k-1).*zstp(1)-zsize./2;
            dist=(dx.*dx+dy.*dy+dz.*dz)^0.5;
            % Density 
            rho(i,j,k)=0;
            if (dist<r)
                rho(i,j,k)=rhoplanet;
            end
            % Right part of Poisson equation
            R1(i,j,k)= 4.0.*G.*pi.*rho(i,j,k);
        end
    end
end

% Defining boundary condition nodes for all grids
for n=levelnum:-1:1
    % Grid points cycle
    for i=1:1:ynum(n);
        for j=1:1:xnum(n);
            for k=1:1:xnum(n);
                % Check distance of (i,j) node from the grid center
                dx=(j-1).*xstp(n)-xsize./2;
                dy=(i-1).*ystp(n)-ysize./2;
                dz=(k-1).*zstp(n)-zsize./2;
                dist=(dx.*dx+dy.*dy+dz.*dz)^0.5;
                % Action depends on the level of resolution
                switch n
                    case 1
                    bon1(i,j,k)=0;
                    if (dist<gradius)
                        bon1(i,j,k)=1;
                    end
                    case 2
                    bon2(i,j,k)=0;
                    if (dist<gradius)
                        bon2(i,j,k)=1;
                    end
                    case 3
                    bon3(i,j,k)=0;
                    if (dist<gradius)
                        bon3(i,j,k)=1;
                    end
                    case 4
                    bon4(i,j,k)=0;
                    if (dist<gradius)
                        bon4(i,j,k)=1;
                    end
                    case 5
                    bon5(i,j,k)=0;
                    if (dist<gradius)
                        bon5(i,j,k)=1;
                    end
                    case 6
                    bon6(i,j,k)=0;
                    if (dist<gradius)
                        bon6(i,j,k)=1;
                    end
                    case 7
                    bon7(i,j,k)=0;
                    if (dist<gradius)
                        bon7(i,j,k)=1;
                    end
                    case 8
                    bon8(i,j,k)=0;
                    if (dist<gradius)
                        bon8(i,j,k)=1;
                    end
                end
            end
        end
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
    for n=1:1:levelnum;
        % Action depends on the level of resolution
        switch n
            
            % Leve 1 (principal grid)
            case 1
            % Smoothing operation: solving of Poisson equation on nodes
            % and computing residuals
            % by calling function Poisson_smoother()
            [phi1,residual1]=Poisson3D_smoother_planet(iternum(n),krelax,xnum(n),ynum(n),znum(n),xstp(n),ystp(n),zstp(n),R1,phi1,bon1,gradius);
            % Restriction operation: 
            % Interpolating residuals to coarcer level (k+1) 
            % to produce right parts for this coarcer level
            if (levelnum>n)
                [R2]=Poisson3D_restriction_planet(n,xnum,ynum,znum,xstp,ystp,zstp,residual1,bon1,bon2);
            end
        
            % Level 2
            case 2
            %Making initial approximation for phi corrections (zeros)
            phi2=zeros(ynum(n),xnum(n),znum(n));
            % Smoothing operation: 
            [phi2,residual2]=Poisson3D_smoother_planet(iternum(n),krelax,xnum(n),ynum(n),znum(n),xstp(n),ystp(n),zstp(n),R2,phi2,bon2,gradius);
            % Restriction operation: 
            if (levelnum>n)
                [R3]=Poisson3D_restriction_planet(n,xnum,ynum,znum,xstp,ystp,zstp,residual2,bon2,bon3);
            end
        
            % Level 3
            case 3
            %Making initial approximation for phi corrections (zeros)
            phi3=zeros(ynum(n),xnum(n),znum(n));
            % Smoothing operation: 
            [phi3,residual3]=Poisson3D_smoother_planet(iternum(n),krelax,xnum(n),ynum(n),znum(n),xstp(n),ystp(n),zstp(n),R3,phi3,bon3,gradius);
            % Restriction operation: 
            if (levelnum>n)
                [R4]=Poisson3D_restriction_planet(n,xnum,ynum,znum,xstp,ystp,zstp,residual3,bon3,bon4);
            end
            
            % Level 4
            case 4
            %Making initial approximation for phi corrections (zeros)
            phi4=zeros(ynum(n),xnum(n),znum(n));
            % Smoothing operation: 
            [phi4,residual4]=Poisson3D_smoother_planet(iternum(n),krelax,xnum(n),ynum(n),znum(n),xstp(n),ystp(n),zstp(n),R4,phi4,bon4,gradius);
            % Restriction operation: 
            if (levelnum>n)
                [R5]=Poisson3D_restriction_planet(n,xnum,ynum,znum,xstp,ystp,zstp,residual4,bon4,bon5);
            end
            
             % Level 5
            case 5
            %Making initial approximation for phi corrections (zeros)
            phi5=zeros(ynum(n),xnum(n),znum(n));
            % Smoothing operation: 
            [phi5,residual5]=Poisson3D_smoother_planet(iternum(n),krelax,xnum(n),ynum(n),znum(n),xstp(n),ystp(n),zstp(n),R5,phi5,bon5,gradius);
            % Restriction operation: 
            if (levelnum>n)
                [R6]=Poisson3D_restriction_planet(n,xnum,ynum,znum,xstp,ystp,zstp,residual5,bon5,bon6);
            end
             
            % Level 6
            case 6
            %Making initial approximation for phi corrections (zeros)
            phi6=zeros(ynum(n),xnum(n),znum(n));
            % Smoothing operation: 
            [phi6,residual6]=Poisson3D_smoother_planet(iternum(n),krelax,xnum(n),ynum(n),znum(n),xstp(n),ystp(n),zstp(n),R6,phi6,bon6,gradius);
              % Restriction operation: 
            if (levelnum>n)
                [R7]=Poisson3D_restriction_planet(n,xnum,ynum,znum,xstp,ystp,zstp,residual6,bon6,bon7);
            end
             
            % Level 7
            case 7
            %Making initial approximation for phi corrections (zeros)
            phi7=zeros(ynum(n),xnum(n),znum(n));
            % Smoothing operation: 
            [phi7,residual7]=Poisson3D_smoother_planet(iternum(n),krelax,xnum(n),ynum(n),znum(n),xstp(n),ystp(n),zstp(n),R7,phi7,bon7,gradius);
              % Restriction operation: 
            if (levelnum>n)
                [R8]=Poisson3D_restriction_planet(n,xnum,ynum,znum,xstp,ystp,zstp,residual7,bon7,bon8);
            end
             
            % Level 8
            case 8
            %Making initial approximation for phi corrections (zeros)
            phi8=zeros(ynum(n),xnum(n),znum(n));
            % Smoothing operation: 
            [phi8,residual8]=Poisson3D_smoother_planet(iternum(n),krelax,xnum(n),ynum(n),znum(n),xstp(n),ystp(n),zstp(n),R8,phi8,bon8,gradius);
            % No Restriction operation on this last level 

        end
    end
   
    % Smoothing+prolongation cycle
    for n=levelnum:-1:1;
        % Action depends on the level of resolution
        switch n

            % Leve 1 (principal grid)
            case 1
            % Smoothing operation: solving of Poisson equation on nodes
            % and computing residuals
            % by calling function Poisson_smoother()
            [phi1,residual1]=Poisson3D_smoother_planet(iternum(n),krelax,xnum(n),ynum(n),znum(n),xstp(n),ystp(n),zstp(n),R1,phi1,bon1,gradius);
            % No prolongation operation on finest level
        
            % Level 2
            case 2
            % Smoothing operation: 
            [phi2,residual2]=Poisson3D_smoother_planet(iternum(n),krelax,xnum(n),ynum(n),znum(n),xstp(n),ystp(n),zstp(n),R2,phi2,bon2,gradius);
            % Prolongation operation: 
            % Interpolating phi corrections to finer level (k+1) 
            % and computing solution update for this finer level
            [dphi1]=Poisson3D_prolongation_planet(n,xnum,ynum,znum,xstp,ystp,zstp,phi2,bon1,bon2);
            phi1=phi1+dphi1;
        
            % Level 3
            case 3
            % Smoothing operation: 
            [phi3,residual3]=Poisson3D_smoother_planet(iternum(n),krelax,xnum(n),ynum(n),znum(n),xstp(n),ystp(n),zstp(n),R3,phi3,bon3,gradius);
            % Prolongation operation: 
            [dphi2]=Poisson3D_prolongation_planet(n,xnum,ynum,znum,xstp,ystp,zstp,phi3,bon2,bon3);
            phi2=phi2+dphi2;
            
            % Level 4
            case 4
            % Smoothing operation: 
            [phi4,residual4]=Poisson3D_smoother_planet(iternum(n),krelax,xnum(n),ynum(n),znum(n),xstp(n),ystp(n),zstp(n),R4,phi4,bon4,gradius);
            % Prolongation operation: 
            [dphi3]=Poisson3D_prolongation_planet(n,xnum,ynum,znum,xstp,ystp,zstp,phi4,bon3,bon4);
            phi3=phi3+dphi3;
            
             % Level 5
            case 5
            % Smoothing operation: 
            [phi5,residual5]=Poisson3D_smoother_planet(iternum(n),krelax,xnum(n),ynum(n),znum(n),xstp(n),ystp(n),zstp(n),R5,phi5,bon5,gradius);
            % Prolongation operation: 
            [dphi4]=Poisson3D_prolongation_planet(n,xnum,ynum,znum,xstp,ystp,zstp,phi5,bon4,bon5);
            phi4=phi4+dphi4;
            
            % Level 6
            case 6
            % Smoothing operation: 
            [phi6,residual6]=Poisson3D_smoother_planet(iternum(n),krelax,xnum(n),ynum(n),znum(n),xstp(n),ystp(n),zstp(n),R6,phi6,bon6,gradius);
            % Prolongation operation: 
            [dphi5]=Poisson3D_prolongation_planet(n,xnum,ynum,znum,xstp,ystp,zstp,phi6,bon5,bon6);
            phi5=phi5+dphi5;
            
            % Level 7
            case 7
            % Smoothing operation: 
            [phi7,residual7]=Poisson3D_smoother_planet(iternum(n),krelax,xnum(n),ynum(n),znum(n),xstp(n),ystp(n),zstp(n),R7,phi7,bon7,gradius);
            % Prolongation operation: 
            [dphi6]=Poisson3D_prolongation_planet(n,xnum,ynum,znum,xstp,ystp,zstp,phi7,bon6,bon7);
            phi6=phi6+dphi6;
            
            % Level 8
            case 8
            % Smoothing operation: 
            [phi8,residual8]=Poisson3D_smoother_planet(iternum(n),krelax,xnum(n),ynum(n),znum(n),xstp(n),ystp(n),zstp(n),R8,phi8,bon8,gradius);
            % Prolongation operation: 
            [dphi7]=Poisson3D_prolongation_planet(n,xnum,ynum,znum,xstp,ystp,zstp,phi8,bon7,bon8);
            phi7=phi7+dphi7;
        end
    end 
 
    % Plotting Residuals for selected stages
    % Normalizing residual
    kfr=4.0.*pi.*G*rhoplanet;
    residual0=residual1/kfr;
    % Defining figure
    figure(1);
    % Plotting Potential
    subplot(2,2,1);
    surf(phi1(:,:,(znum(1)-1)/2));
    %colormap(Gray);
    shading interp;
    light;
    lighting phong;
    axis tight;
    zlabel('Gravity potential');
    title(['Solution of Poisson equation, V-cycle = ',num2str(niter)]);
    % Plotting residual
    subplot(2,2,2);
    surf(residual0(:,:,(znum(1)-1)/2));
    %colormap(Gray);
    shading interp;
    light;
    lighting phong;
    axis tight;
    zlabel('residual');
    title(['Solution of Poisson equation, V-cycle = ',num2str(niter)]);
    
    % Computing mean square residuals
    resphi00(niter)=0;
    for i=1:1:ynum(1)
        for j=1:1:xnum(1)
            for k=1:1:znum(1)
                resphi00(niter)=resphi00(niter)+residual0(i,j,k)^2;
            end
        end
    end
    % Compute log of new  mean square residuals
    resphi00(niter)=log10((resphi00(niter)/(ynum(1)*xnum(1)*znum(1)))^0.5);
    % Plotting Mean residuals
    subplot(2,2,3)
    plot(resphi00, 'k');
    xlabel('V-cycles')
    ylabel(['log(Residuals) iteraton ',num2str(niter)])
       

    pause(3);

end

