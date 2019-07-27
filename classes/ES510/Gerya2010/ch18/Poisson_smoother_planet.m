% Function Poisson_smoother()
% This function makes specified number (iternum) of Gauss Seidel iteraions 
% using relaxation koefficient (krelax)
% for Poisson equation defined on 2D grid
% with specified resolution (xnum, ynum) and gridsteps (xstp, ystp)
% given distribution of right parts for Poisson equations (R) on the grid 
% and given approximation of gravity potential (phi) on the grid 
%
% Function return new approximation for the gravity potential (phi)
% and distribution of residuals (residual)
function[phi,residual]=Poisson_smoother_planet(iternum,krelax,xnum,ynum,xstp,ystp,R,phi,bon,gradius)

% Poisson equation koefficients
xkf=1/xstp^2;
ykf=1/ystp^2;
xykf=-2/xstp^2-2/ystp^2;
% Planet center coordinates
xcenter=(xnum-1)*xstp/2;
ycenter=(ynum-1)*ystp/2;
% Square of radius of boundary circle
gradius2=gradius*gradius;

% Gauss-Seidel iteration circle
for niter=1:1:iternum;
    % Solving of Poisson equation on nodes
    for i=1:1:ynum;
        for j=1:1:xnum;
            % Boundrary conditions 
            if (bon(i,j)==0 || i==1 || i==ynum || j==1 || j==xnum);
                % Boundary condition for gravity potential
                phi(i,j)=0;
            % Solving Poisson equation
            else
                % Initial coefficient value for the central node
                kfcur=xykf;
                rescur=R(i,j);
                % Check for surrounding ghost Nodes outside of boundary condition circle
                % Ghost node to the left
                if (bon(i,j-1)==0)
                    % Compute horizontal distance to the boundary circle
                    dx=abs((j-1)*xstp-xcenter);
                    dy=abs((i-1)*ystp-ycenter);
                    dxr=(gradius2-dy*dy)^0.5-dx;
                    % Compute koefficient for boundary condition equation
                    % phi(i,j-1)=-phi(i,j)*(xstp-dxr)/dxr
                    % satisfying phi=0 on the boundary circle
                    kfcur=kfcur-xkf*(xstp-dxr)/dxr;
                % Normal node to the left
                else
                    rescur=rescur-xkf*phi(i,j-1);
                end
                % Ghost node to the right
                if (bon(i,j+1)==0)
                    % Compute horizontal distance to the boundary circle
                    dx=abs((j-1)*xstp-xcenter);
                    dy=abs((i-1)*ystp-ycenter);
                    dxr=(gradius2-dy*dy)^0.5-dx;
                    % Compute koefficient for boundary condition equation
                    % phi(i,j+1)=-phi(i,j)*(xstp-dxr)/dxr
                    % satisfying phi=0 on the boundary circle
                    kfcur=kfcur-xkf*(xstp-dxr)/dxr;
                % Normal node to the right
                else
                    rescur=rescur-xkf*phi(i,j+1);
                end
                % Ghost node to the top
                if (bon(i-1,j)==0)
                    % Compute vertical distance to the boundary circle
                    dx=abs((j-1)*xstp-xcenter);
                    dy=abs((i-1)*ystp-ycenter);
                    dyr=(gradius2-dx*dx)^0.5-dy;
                    % Compute koefficient for boundary condition equation
                    % phi(i-1,j)=-phi(i,j)*(xstp-dyr)/dyr
                    % satisfying phi=0 on the boundary circle
                    kfcur=kfcur-ykf*(ystp-dyr)/dyr;
                % Normal node to the top
                else
                    rescur=rescur-ykf*phi(i-1,j);
                end
                % Ghost node to the bottom
                if (bon(i+1,j)==0)
                    % Compute vertical distance to the boundary circle
                    dx=abs((j-1)*xstp-xcenter);
                    dy=abs((i-1)*ystp-ycenter);
                    dyr=(gradius2-dx*dx)^0.5-dy;
                    % Compute koefficient for boundary condition equation
                    % phi(i+1,j)=-phi(i,j)*(xstp-dyr)/dyr
                    % satisfying phi=0 on the boundary circle
                    kfcur=kfcur-ykf*(ystp-dyr)/dyr;
                % Normal node to the top
                else
                    rescur=rescur-ykf*phi(i+1,j);
                end
                % Computing current residual
                rescur=rescur-kfcur*phi(i,j);
                % Updating solution
                phi(i,j)=phi(i,j)+rescur/kfcur*krelax;
            end
        end            
    end
    % End of Solving Poisson equation on nodes
end
% Computing final state of residuals
for i=1:1:ynum;
    for j=1:1:xnum;
        % Residuals for boundary conditions 
        if (bon(i,j)==0 || i==1 || i==ynum || j==1 || j==xnum);
            residual(i,j)=0;
        else
            % Residuals for Poisson equation
            % Initial coefficient value for the central node
            kfcur=xykf;
            rescur=R(i,j);
            % Check for surrounding ghost Nodes outside of boundary condition circle
            % Ghost node to the left
            if (bon(i,j-1)==0)
                % Compute horizontal distance to the boundary circle
                dx=abs((j-1)*xstp-xcenter);
                dy=abs((i-1)*ystp-ycenter);
                dxr=(gradius2-dy*dy)^0.5-dx;
                % Compute koefficient for boundary condition equation
                % phi(i,j-1)=-phi(i,j)*(xstp-dxr)/dxr
                % satisfying phi=0 on the boundary circle
                kfcur=kfcur-xkf*(xstp-dxr)/dxr;
            % Normal node to the left
            else
                rescur=rescur-xkf*phi(i,j-1);
            end
            % Ghost node to the right
            if (bon(i,j+1)==0)
                % Compute horizontal distance to the boundary circle
                dx=abs((j-1)*xstp-xcenter);
                dy=abs((i-1)*ystp-ycenter);
                dxr=(gradius2-dy*dy)^0.5-dx;
                % Compute koefficient for boundary condition equation
                % phi(i,j+1)=-phi(i,j)*(xstp-dxr)/dxr
                % satisfying phi=0 on the boundary circle
                kfcur=kfcur-xkf*(xstp-dxr)/dxr;
            % Normal node to the right
            else
                rescur=rescur-xkf*phi(i,j+1);
            end
            % Ghost node to the top
            if (bon(i-1,j)==0)
                % Compute vertical distance to the boundary circle
                dx=abs((j-1)*xstp-xcenter);
                dy=abs((i-1)*ystp-ycenter);
                dyr=(gradius2-dx*dx)^0.5-dy;
                % Compute koefficient for boundary condition equation
                % phi(i-1,j)=-phi(i,j)*(xstp-dyr)/dyr
                % satisfying phi=0 on the boundary circle
                kfcur=kfcur-ykf*(ystp-dyr)/dyr;
            % Normal node to the top
            else
                rescur=rescur-ykf*phi(i-1,j);
            end
            % Ghost node to the bottom
            if (bon(i+1,j)==0)
                % Compute vertical distance to the boundary circle
                dx=abs((j-1)*xstp-xcenter);
                dy=abs((i-1)*ystp-ycenter);
                dyr=(gradius2-dx*dx)^0.5-dy;
                % Compute koefficient for boundary condition equation
                % phi(i+1,j)=-phi(i,j)*(xstp-dyr)/dyr
                % satisfying phi=0 on the boundary circle
                kfcur=kfcur-ykf*(ystp-dyr)/dyr;
            % Normal node to the top
            else
                rescur=rescur-ykf*phi(i+1,j);
            end
            % Computing current residual
            residual(i,j)=rescur-kfcur*phi(i,j);
        end
    end            
end

