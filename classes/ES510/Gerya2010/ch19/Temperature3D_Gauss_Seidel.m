% Solution of Temperature equation in 3D with Gauss-Seidel iteration
% by using external function Temperature3D_smoother()
% Temperature and density distribution in the model  
% corresponds to a hot dense less conductive block in the cold lower density medium

% Clearing all variables and arrays
clear all;
% Clearing figures
clf;

% Model parameters
% Block density
rhoblock=3100.0;
% Block thermal conductivity W*m/K
ktblock=1;
% Heat capacity J/K/kg
cpblock=1500;
% Block temperature, K
tkblock=1500;

% Medium density
rhomedium=3000.0;
% Block thermal conductivity W*m/K
ktmedium=3;
% Heat capacity J/K/kg
cpmedium=1000;
% Medium temperature, K
tkmedium=1000;

% Model size, m
xsize=100000.0;
ysize=100000.0;
zsize=100000.0;

% Defining resolutions on all levels
xnum=51;
ynum=51;
znum=51;

% Grid steps
xstp=xsize./(xnum-1);
ystp=ysize./(ynum-1);
zstp=zsize./(znum-1);

% Number of timesteps
ntimesteps=100;
% Timestep, yr
tstep=1e+5;
% Timestep, s
tstep=tstep*(365.25*24*3600);

% Number of smoothing iterations before checking residuals
iternum=5;
% Required accuracy of solving
taccuracy=1e-09;
% Minimal number of iterations
nitermin=15;

% Relaxation coefficients for Gauss-Seidel iterations
krelax=1.25;


% FINEST (PRINCIPAL) GRID
% Initial teperature structure tk()
% Defining distribution of density rho(), thermal conductivity kt(), heat capacity cp(), 
% Grid points cycle
for i=1:1:ynum;
    for j=1:1:xnum;
        for k=1:1:znum;

            % Relative distances for (i,j,k) basic node inside the grid
            dx=(j-1)*xstp(1)/xsize;
            dy=(i-1)*ystp(1)/ysize;
            dz=(k-1)*zstp(1)/zsize;
            % T, rho, cp, kt Structure (block in a medium)
            tk1(i,j,k)=tkmedium;
            rho(i,j,k)=rhomedium;
            kt(i,j,k)=ktmedium;
            cp(i,j,k)=cpmedium;
            if(dx>=0.4 && dx<=0.6 && dy>=0.4 && dy<=0.6 && dz>=0.4 && dz<=0.6)
                tk1(i,j,k)=tkblock;
                rho(i,j,k)=rhoblock;
                kt(i,j,k)=ktblock;
                cp(i,j,k)=cpblock;
            end
            
        end
    end
end

% Reset number of iterations 
niter=1;
%Reset time
timesum=0;

% Time cycle
for nnn=1:1:ntimesteps;
    
    % Reset current accuracy
    taccur=taccuracy*1e+20;
    % Save current temperature
    tk0=tk1;
    % Reset counter of iterations per timestep
    nitercur=0;

    % Solve for the new temperature
    while (taccur>taccuracy && nitercur<nitermin)
        % Smoothing operation: solving of Temperature equation on nodes
        % and computing residuals
        % by calling function Temperature3D_smoother()
        [tk1,residual1]=Temperature3D_smoother(iternum,krelax,xnum,ynum,znum,xstp,ystp,zstp,tk0,tk1,rho,cp,kt,tstep);
 
        % Plotting Residuals for selected stages
        % Normalizing residual
        kfr=tkmedium*cpmedium*rhomedium/tstep;
        residual0=residual1/kfr;
        
        % Defining figure
        figure(1);
        
        % Plotting New Temperature
        subplot(2,2,1);
        surf(tk1(:,:,(znum-1)/2));
        %colormap(Gray);
        shading interp;
        light;
        lighting phong;
        axis tight;
        zlabel(['T K, time Myr ',num2str(timesum/(365.25*24*3600*1e+6))]);
%       title(['Solution of Temperature equation, Gauss-Seidel cycle = ',num2str(niter)]);

        % Plotting Tnew-Told
        subplot(2,2,2);
        dtk1=tk1-tk0;
        surf(dtk1(:,:,(znum-1)/2));
        %colormap(Gray);
        shading interp;
        light;
        lighting phong;
        axis tight;
        zlabel('dT K');
%       title(['Temperature equation, Gauss-Seidel cycle = ',num2str(niter)]);
    
        % Plotting residual
        subplot(2,2,3);
        surf(residual0(:,:,(znum-1)/2));
        %colormap(Gray);
        shading interp;
        light;
        lighting phong;
        axis tight;
        zlabel('residual');
%       title(['Solution of Temperature equation, Gauss-Seidel cycle = ',num2str(niter)]);
    
        % Computing mean square residuals
        rest00(niter)=0;
        for i=1:1:ynum
            for j=1:1:xnum
                for k=1:1:znum
                    rest00(niter)=rest00(niter)+residual0(i,j,k)^2;
                end
            end
        end
        % Compute log of new  mean square residuals
        rest00(niter)=log10((rest00(niter)/(ynum*xnum*znum))^0.5);
        % Save new accuracy
        taccur=10^rest00(niter);
 
        % Plotting Mean residuals
        subplot(2,2,4)
        plot(rest00, 'k');
        xlabel('Gauss-Seidel-cycles')
        ylabel(['log(Residuals) iteraton ',num2str(niter)])
        
        % Add iterations for the timestep counter
        nitercur=nitercur+1;
        % Add total iterations counter
        niter=niter+1;

    
        pause(1);
    end
    % Add time
    timesum=timesum+tstep;

end

