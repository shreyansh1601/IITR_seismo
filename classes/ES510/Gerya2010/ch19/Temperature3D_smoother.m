% Function Temperature3D_smoother()
% This function makes specified number (iternum) of Gauss Seidel iteraions 
% using relaxation koefficient (krelax)
% for Poisson equation defined on 2D grid
% with specified resolution (xnum, ynum) and gridsteps (xstp, ystp)
% given distribution of right parts for Poisson equations (R) on the grid 
% and given approximation of gravity potential (phi) on the grid 
%
% Function return new approximation for the gravity potential (phi)
% and distribution of residuals (residual)
function[tk1,residual]=Temperature3D_smoother(iternum,krelax,xnum,ynum,znum,xstp,ystp,zstp,tk0,tk1,rho,cp,kt,tstep)

% Poisson equation koefficients
xkf=1/xstp^2/2;
ykf=1/ystp^2/2;
zkf=1/zstp^2/2;

% Gauss-Seidel iteration circle
for niter=1:1:iternum;
    % Solving of Poisson equation on nodes
    for i=1:1:ynum;
        for j=1:1:xnum;
            for k=1:1:znum;
                % Boundary conditions 
                if (i==1 || i==ynum || j==1 || j==xnum || k==1 || k==znum)
                     % Top boundary: symmetry
                    if (i==1)
                        tk1(i,j,k)=tk1(i+1,j,k);
                    end
                    % Bottom boundary: symmetry
                    if (i==ynum)
                        tk1(i,j,k)=tk1(i-1,j,k);
                    end
                   % Left boundary: symmetry
                    if (j==1)
                        tk1(i,j,k)=tk1(i,j+1,k);
                    end
                    % Right boundary: symmetry
                    if (j==xnum)
                        tk1(i,j,k)=tk1(i,j-1,k);
                    end
                    % Back boundary: symmetry
                    if (k==1)
                        tk1(i,j,k)=tk1(i,j,k+1);
                    end
                    % Front boundary: symmetry
                    if (k==znum)
                        tk1(i,j,k)=tk1(i,j,k-1);
                    end
                % Solving temperature equation
                else
                    % Residuals for Temperature equation
                    % dT/dt
                    rescur=rho(i,j,k)*cp(i,j,k)*(tk0(i,j,k)-tk1(i,j,k))/tstep;
                    % -dqx/dx
                    rescur=rescur+xkf*((kt(i,j,k)+kt(i,j+1,k))*(tk1(i,j+1,k)-tk1(i,j,k))-(kt(i,j-1,k)+kt(i,j,k))*(tk1(i,j,k)-tk1(i,j-1,k)));
                    % -dqy/dy
                    rescur=rescur+ykf*((kt(i,j,k)+kt(i+1,j,k))*(tk1(i+1,j,k)-tk1(i,j,k))-(kt(i-1,j,k)+kt(i,j,k))*(tk1(i,j,k)-tk1(i-1,j,k)));
                    % -dqz/dz
                    rescur=rescur+zkf*((kt(i,j,k)+kt(i,j,k+1))*(tk1(i,j,k+1)-tk1(i,j,k))-(kt(i,j,k-1)+kt(i,j,k))*(tk1(i,j,k)-tk1(i,j,k-1)));
                    % Coefficient value for the central node
                    kfcur=rho(i,j,k)*cp(i,j,k)/tstep+xkf*(kt(i,j,k)+kt(i,j+1,k)+kt(i,j-1,k)+kt(i,j,k))+ykf*(kt(i,j,k)+kt(i+1,j,k)+kt(i-1,j,k)+kt(i,j,k))+zkf*(kt(i,j,k)+kt(i,j,k+1)+kt(i,j,k-1)+kt(i,j,k));
                    % Updating solution
                    tk1(i,j,k)=tk1(i,j,k)+rescur/kfcur*krelax;
                end
            end
        end            
    end
    % End of Solving Temperature equation on nodes
end

% Computing final state of residuals
for i=1:1:ynum;
    for j=1:1:xnum;
        for k=1:1:znum;
            % Residuals for boundary conditions 
            if (i==1 || i==ynum || j==1 || j==xnum || k==1 || k==znum)
                residual(i,j,k)=0;
            % Residual for temperature equation
            else
                % Residuals for Temperature equation
                % dT/dt
                rescur=rho(i,j,k)*cp(i,j,k)*(tk0(i,j,k)-tk1(i,j,k))/tstep;
                % -dqx/dx
                rescur=rescur+xkf*((kt(i,j,k)+kt(i,j+1,k))*(tk1(i,j+1,k)-tk1(i,j,k))-(kt(i,j-1,k)+kt(i,j,k))*(tk1(i,j,k)-tk1(i,j-1,k)));
                % -dqy/dy
                rescur=rescur+ykf*((kt(i,j,k)+kt(i+1,j,k))*(tk1(i+1,j,k)-tk1(i,j,k))-(kt(i-1,j,k)+kt(i,j,k))*(tk1(i,j,k)-tk1(i-1,j,k)));
                % -dqz/dz
                rescur=rescur+zkf*((kt(i,j,k)+kt(i,j,k+1))*(tk1(i,j,k+1)-tk1(i,j,k))-(kt(i,j,k-1)+kt(i,j,k))*(tk1(i,j,k)-tk1(i,j,k-1)));
                residual(i,j,k)=rescur;
            end
        end
    end            
end

