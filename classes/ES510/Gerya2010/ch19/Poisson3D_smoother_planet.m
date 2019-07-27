% Function Poisson3D_smoother_planet()
% This function makes specified number (iternum) of Gauss Seidel iteraions 
% using relaxation koefficient (krelax)
% for Poisson equation defined on 2D grid
% with specified resolution (xnum, ynum) and gridsteps (xstp, ystp)
% given distribution of right parts for Poisson equations (R) on the grid 
% and given approximation of gravity potential (phi) on the grid 
%
% Function return new approximation for the gravity potential (phi)
% and distribution of residuals (residual)
function[phi,residual]=Poisson3D_smoother_planet(iternum,krelax,xnum,ynum,znum,xstp,ystp,zstp,R,phi,bon,gradius)

% Poisson equation koefficients
xkf=1/xstp^2;
ykf=1/ystp^2;
zkf=1/zstp^2;
xyzkf=-2/xstp^2-2/ystp^2-2/zstp^2;
% Planet center coordinates
xcenter=(xnum-1)*xstp/2;
ycenter=(ynum-1)*ystp/2;
zcenter=(znum-1)*zstp/2;
% Square of radius of boundary circle
gradius2=gradius*gradius;

% Gauss-Seidel iteration circle
for niter=1:1:iternum;
    % Solving of Poisson equation on nodes
    for i=1:1:ynum;
        for j=1:1:xnum;
            for k=1:1:znum;
            % Residuals for boundary conditions 
            if (bon(i,j,k)==0 || i==1 || i==ynum || j==1 || j==xnum || k==1 || k==znum);
                phi(i,j,k)=0;
            else
                % Residuals for Poisson equation
                % Initial coefficient value for the central node
                kfcur=xyzkf;
                rescur=R(i,j,k);
                % Check for surrounding ghost Nodes outside of boundary condition circle
                % Ghost node to the left
                if (bon(i,j-1,k)==0)
                    % Compute horizontal distance to the boundary circle
                    dx=abs((j-1)*xstp-xcenter);
                    dy=abs((i-1)*ystp-ycenter);
                    dz=abs((k-1)*zstp-zcenter);
                    dxr=(gradius2-dy*dy-dz*dz)^0.5-dx;
                    % Compute koefficient for boundary condition equation
                    % phi(i,j-1,k)=-phi(i,j,k)*(xstp-dxr)/dxr
                    % satisfying phi=0 on the boundary circle
                    kfcur=kfcur-xkf*(xstp-dxr)/dxr;
                % Normal node to the left
                else
                    rescur=rescur-xkf*phi(i,j-1,k);
                end
                % Ghost node to the right
                if (bon(i,j+1,k)==0)
                    % Compute horizontal distance to the boundary circle
                    dx=abs((j-1)*xstp-xcenter);
                    dy=abs((i-1)*ystp-ycenter);
                    dz=abs((k-1)*zstp-zcenter);
                    dxr=(gradius2-dy*dy-dz*dz)^0.5-dx;
                    % Compute koefficient for boundary condition equation
                    % phi(i,j+1,k)=-phi(i,j,k)*(xstp-dxr)/dxr
                    % satisfying phi=0 on the boundary circle
                    kfcur=kfcur-xkf*(xstp-dxr)/dxr;
                % Normal node to the right
                else
                    rescur=rescur-xkf*phi(i,j+1,k);
                end
                % Ghost node to the top
                if (bon(i-1,j,k)==0)
                    % Compute vertical distance to the boundary circle
                    dx=abs((j-1)*xstp-xcenter);
                    dy=abs((i-1)*ystp-ycenter);
                    dz=abs((k-1)*zstp-zcenter);
                    dyr=(gradius2-dx*dx-dz*dz)^0.5-dy;
                    % Compute koefficient for boundary condition equation
                    % phi(i-1,j,k)=-phi(i,j,k)*(ystp-dyr)/dyr
                    % satisfying phi=0 on the boundary circle
                    kfcur=kfcur-ykf*(ystp-dyr)/dyr;
                % Normal node to the top
                else
                    rescur=rescur-ykf*phi(i-1,j,k);
                end
                % Ghost node to the bottom
                if (bon(i+1,j,k)==0)
                    % Compute vertical distance to the boundary circle
                    dx=abs((j-1)*xstp-xcenter);
                    dy=abs((i-1)*ystp-ycenter);
                    dz=abs((k-1)*zstp-zcenter);
                    dyr=(gradius2-dx*dx-dz*dz)^0.5-dy;
                    % Compute koefficient for boundary condition equation
                    % phi(i+1,j,k)=-phi(i,j,k)*(ystp-dyr)/dyr
                    % satisfying phi=0 on the boundary circle
                    kfcur=kfcur-ykf*(ystp-dyr)/dyr;
                % Normal node to the bottom
                else
                    rescur=rescur-ykf*phi(i+1,j,k);
                end
                % Ghost node to the back
                if (bon(i,j,k-1)==0)
                    % Compute vertical distance to the boundary circle
                    dx=abs((j-1)*xstp-xcenter);
                    dy=abs((i-1)*ystp-ycenter);
                    dz=abs((k-1)*zstp-zcenter);
                    dzr=(gradius2-dx*dx-dy*dy)^0.5-dz;
                    % Compute koefficient for boundary condition equation
                    % phi(i,j,k-1)=-phi(i,j,k)*(zstp-dzr)/dzr
                    % satisfying phi=0 on the boundary circle
                    kfcur=kfcur-zkf*(zstp-dzr)/dzr;
                % Normal node to the back
                else
                    rescur=rescur-zkf*phi(i,j,k-1);
                end
                % Ghost node to the bottom
                if (bon(i,j,k+1)==0)
                    % Compute vertical distance to the boundary circle
                    dx=abs((j-1)*xstp-xcenter);
                    dy=abs((i-1)*ystp-ycenter);
                    dz=abs((k-1)*zstp-zcenter);
                    dzr=(gradius2-dx*dx-dy*dy)^0.5-dz;
                    % Compute koefficient for boundary condition equation
                    % phi(i,j,k+1)=-phi(i,j,k)*(zstp-dzr)/dzr
                    % satisfying phi=0 on the boundary circle
                    kfcur=kfcur-zkf*(zstp-dzr)/dzr;
                % Normal node to the front
                else
                    rescur=rescur-zkf*phi(i,j,k+1);
                end
                    % Computing current residual
                    rescur=rescur-kfcur*phi(i,j,k);
                    % Updating solution
                    phi(i,j,k)=phi(i,j,k)+rescur/kfcur*krelax;
                end
            end
        end            
    end
    % End of Solving Poisson equation on nodes
end
% Computing final state of residuals
for i=1:1:ynum;
    for j=1:1:xnum;
        for k=1:1:znum;
            % Residuals for boundary conditions 
            if (bon(i,j,k)==0 || i==1 || i==ynum || j==1 || j==xnum || k==1 || k==znum);
                residual(i,j,k)=0;
            else
                % Residuals for Poisson equation
                % Initial coefficient value for the central node
                kfcur=xyzkf;
                rescur=R(i,j,k);
                % Check for surrounding ghost Nodes outside of boundary condition circle
                % Ghost node to the left
                if (bon(i,j-1,k)==0)
                    % Compute horizontal distance to the boundary circle
                    dx=abs((j-1)*xstp-xcenter);
                    dy=abs((i-1)*ystp-ycenter);
                    dz=abs((k-1)*zstp-zcenter);
                    dxr=(gradius2-dy*dy-dz*dz)^0.5-dx;
                    % Compute koefficient for boundary condition equation
                    % phi(i,j-1,k)=-phi(i,j,k)*(xstp-dxr)/dxr
                    % satisfying phi=0 on the boundary circle
                    kfcur=kfcur-xkf*(xstp-dxr)/dxr;
                % Normal node to the left
                else
                    rescur=rescur-xkf*phi(i,j-1,k);
                end
                % Ghost node to the right
                if (bon(i,j+1,k)==0)
                    % Compute horizontal distance to the boundary circle
                    dx=abs((j-1)*xstp-xcenter);
                    dy=abs((i-1)*ystp-ycenter);
                    dz=abs((k-1)*zstp-zcenter);
                    dxr=(gradius2-dy*dy-dz*dz)^0.5-dx;
                    % Compute koefficient for boundary condition equation
                    % phi(i,j+1,k)=-phi(i,j,k)*(xstp-dxr)/dxr
                    % satisfying phi=0 on the boundary circle
                    kfcur=kfcur-xkf*(xstp-dxr)/dxr;
                % Normal node to the right
                else
                    rescur=rescur-xkf*phi(i,j+1,k);
                end
                % Ghost node to the top
                if (bon(i-1,j,k)==0)
                    % Compute vertical distance to the boundary circle
                    dx=abs((j-1)*xstp-xcenter);
                    dy=abs((i-1)*ystp-ycenter);
                    dz=abs((k-1)*zstp-zcenter);
                    dyr=(gradius2-dx*dx-dz*dz)^0.5-dy;
                    % Compute koefficient for boundary condition equation
                    % phi(i-1,j,k)=-phi(i,j,k)*(ystp-dyr)/dyr
                    % satisfying phi=0 on the boundary circle
                    kfcur=kfcur-ykf*(ystp-dyr)/dyr;
                % Normal node to the top
                else
                    rescur=rescur-ykf*phi(i-1,j,k);
                end
                % Ghost node to the bottom
                if (bon(i+1,j,k)==0)
                    % Compute vertical distance to the boundary circle
                    dx=abs((j-1)*xstp-xcenter);
                    dy=abs((i-1)*ystp-ycenter);
                    dz=abs((k-1)*zstp-zcenter);
                    dyr=(gradius2-dx*dx-dz*dz)^0.5-dy;
                    % Compute koefficient for boundary condition equation
                    % phi(i+1,j,k)=-phi(i,j,k)*(ystp-dyr)/dyr
                    % satisfying phi=0 on the boundary circle
                    kfcur=kfcur-ykf*(ystp-dyr)/dyr;
                % Normal node to the bottom
                else
                    rescur=rescur-ykf*phi(i+1,j,k);
                end
                % Ghost node to the back
                if (bon(i,j,k-1)==0)
                    % Compute vertical distance to the boundary circle
                    dx=abs((j-1)*xstp-xcenter);
                    dy=abs((i-1)*ystp-ycenter);
                    dz=abs((k-1)*zstp-zcenter);
                    dzr=(gradius2-dx*dx-dy*dy)^0.5-dz;
                    % Compute koefficient for boundary condition equation
                    % phi(i,j,k-1)=-phi(i,j,k)*(zstp-dzr)/dzr
                    % satisfying phi=0 on the boundary circle
                    kfcur=kfcur-zkf*(zstp-dzr)/dzr;
                % Normal node to the back
                else
                    rescur=rescur-zkf*phi(i,j,k-1);
                end
                % Ghost node to the bottom
                if (bon(i,j,k+1)==0)
                    % Compute vertical distance to the boundary circle
                    dx=abs((j-1)*xstp-xcenter);
                    dy=abs((i-1)*ystp-ycenter);
                    dz=abs((k-1)*zstp-zcenter);
                    dzr=(gradius2-dx*dx-dy*dy)^0.5-dz;
                    % Compute koefficient for boundary condition equation
                    % phi(i,j,k+1)=-phi(i,j,k)*(zstp-dzr)/dzr
                    % satisfying phi=0 on the boundary circle
                    kfcur=kfcur-zkf*(zstp-dzr)/dzr;
                % Normal node to the front
                else
                    rescur=rescur-zkf*phi(i,j,k+1);
                end
                % Computing current residual
                residual(i,j,k)=rescur-kfcur*phi(i,j,k);
            end
        end
    end            
end

