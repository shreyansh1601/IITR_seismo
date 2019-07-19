% Function Poisson_solver_grid()
% This function formulates and solves  
% Poisson equation defined on 2D irregularly spaced grid
% with specified resolution (xnum, ynum) and grid lines positions (gridx, gridy)
% given radius (rbon) and value (vbon) for boundary condition domain 
% from given planetary centre (xcentre,ycentre) 
% given distribution of right parts for all equations (RP) on the grid 
%
% Function returns solution for new potential (finew)
% and distribution of residuals (resfi)
function[L,R,finew,resfi]=Poisson_solver_planet_grid(L,R,xnum,ynum,gridx,gridy,RP,xcentre,ycentre,rbon,vbon)

% Computing grid steps for FI nodes
xstp=zeros(xnum-1,1);
ystp=zeros(ynum-1,1);
for i=1:1:xnum-1
    xstp(i)=gridx(i+1)-gridx(i);
end
for i=1:1:ynum-1
    ystp(i)=gridy(i+1)-gridy(i);
end


% Solving of Poisson equation on nodes
for i=1:1:ynum
    for j=1:1:xnum
        % Index for G
        itk=(j-1)*ynum+i;
        
        % Boundary conditions
        % Check distance from planetary centre
        dx=gridx(j)-xcentre;
        dy=gridy(i)-ycentre;
        dr=(dx^2+dy^2)^0.5;
        if (j==1 || j==xnum || i==1 || i==ynum || dr>=rbon)
            % Constant given value
            % Right part
            R(itk,1)=vbon;
            % Left part: 1*fi(i,j)
            L(itk,itk)=1;
        % Poisson equation d2FI/dx2+d2FI/dy2=KP*4*PI*GAMMA*RHO
        else
            % Poisson equation stensil
            %
            %     +-------fi(i-1,j)-------+-
            %     |           |           |
            %     |           |           |  
            %     |           |           |
            %     |           |           |
            % fi(i,j-1)----fi(i,j)-----fi(i,j+1)
            %     |           |           |  
            %     |           |           |
            %     |           |           |
            %     +-------fi(i+1,j)-------+-
            %
            % Right Part
            R(itk,1)=RP(i,j);
            % Computing coefficients for the left part
            L(itk,itk)=0; % Initial coefficient for fi(i,j)
            % dFI/dx(j-1/2)
            % Check distance from fi(i,j-1)to the planetary centre
             % Check distance from planetary centre
            dx1=gridx(j-1)-xcentre;
            dr1=(dx1^2+dy^2)^0.5;
            % fi(i,j-1) node is in the boundary condition domain
            if (dr1>=rbon)
                % Compute distance to the external boundary
                dx1=(rbon^2-dy^2)^0.5-abs(dx);
                % Add dFI/dx derivative to the left part
                L(itk,itk)=L(itk,itk)-1/dx1/((xstp(j-1)+xstp(j))/2);
                % Add dFI/dx derivative to the right part
                R(itk,1)=R(itk,1)-vbon/dx1/((xstp(j-1)+xstp(j))/2);
            else
            % gp(i,j-1) node is not in the boundary condition domain
                % Add dFI/dx derivative to the left part
                L(itk,itk)=L(itk,itk)-1/xstp(j-1)/((xstp(j-1)+xstp(j))/2); % gp(i,j) node
                L(itk,itk-ynum)= 1/xstp(j-1)/((xstp(j-1)+xstp(j))/2); % gp(i,j-1) node
            end
            % -dFI/dx(j+1/2)
            % Check distance from fi(i,j+1)to the planetary centre
             % Check distance from planetary centre
            dx1=gridx(j+1)-xcentre;
            dr1=(dx1^2+dy^2)^0.5;
            % fi(i,j-1) node is in the boundary condition domain
            if (dr1>=rbon)
                % Compute distance to the external boundary
                dx1=(rbon^2-dy^2)^0.5-abs(dx);
                % Add dFI/dx derivative to the left part
                L(itk,itk)=L(itk,itk)-1/dx1/((xstp(j-1)+xstp(j))/2);
                % Add dFI/dx derivative to the right part
                R(itk,1)=R(itk,1)-vbon/dx1/((xstp(j-1)+xstp(j))/2);
            else
            % fi(i,j+1) node is not in the boundary condition domain
                % Add dFI/dx derivative to the left part
                L(itk,itk)=L(itk,itk)-1/xstp(j)/((xstp(j-1)+xstp(j))/2); % gp(i,j) node
                L(itk,itk+ynum)= 1/xstp(j)/((xstp(j-1)+xstp(j))/2); % gp(i,j-1) node
            end
            % -dFI/dy(i-1/2)
            % Check distance from fi(i-1,j)to the planetary centre
             % Check distance from planetary centre
            dy1=gridy(i-1)-ycentre;
            dr1=(dx^2+dy1^2)^0.5;
            % fi(i-1,j) node is in the boundary condition domain
            if (dr1>=rbon)
                % Compute distance to the external boundary
                dy1=(rbon^2-dx^2)^0.5-abs(dy);
                % Add dFI/dx derivative to the left part
                L(itk,itk)=L(itk,itk)-1/dy1/((ystp(i-1)+ystp(i))/2);
                % Add dFI/dx derivative to the right part
                R(itk,1)=R(itk,1)-vbon/dy1/((ystp(i-1)+ystp(i))/2);
            else
            % fi(i-1,j) node is not in the boundary condition domain
                % Add dFI/dx derivative to the left part
                L(itk,itk)=L(itk,itk)-1/ystp(i-1)/((ystp(i-1)+ystp(i))/2); % gp(i,j) node
                L(itk,itk-1)= 1/ystp(i-1)/((ystp(i-1)+ystp(i))/2); % gp(i-1,j) node
            end
            % dFI/dy(i+1/2)
            % Check distance from fi(i+1,j)to the planetary centre
             % Check distance from planetary centre
            dy1=gridy(i+1)-ycentre;
            dr1=(dx^2+dy1^2)^0.5;
            % fi(i+1,j) node is in the boundary condition domain
            if (dr1>=rbon)
                % Compute distance to the external boundary
                dy1=(rbon^2-dx^2)^0.5-abs(dy);
                % Add dFI/dx derivative to the left part
                L(itk,itk)=L(itk,itk)-1/dy1/((ystp(i-1)+ystp(i))/2);
                % Add dFI/dx derivative to the right part
                R(itk,1)=R(itk,1)-vbon/dy1/((ystp(i-1)+ystp(i))/2);
            else
            % fi(i+1,j) node is not in the boundary condition domain
                % Add dFI/dx derivative to the left part
                L(itk,itk)=L(itk,itk)-1/ystp(i)/((ystp(i-1)+ystp(i))/2); % gp(i,j) node
                L(itk,itk+1)= 1/ystp(i)/((ystp(i-1)+ystp(i))/2); % gp(i-1,j) node
            end
        end
    end            
end


% Solve matrix
S=L\R;
finew=zeros(ynum,xnum);
% Reload solution
for i=1:1:ynum
    for j=1:1:xnum
        % Index for T
        itk=(j-1)*ynum+i;
        % Reload T
        finew(i,j)=S(itk);
    end
end


% Computing residuals
resfi=zeros(ynum,xnum);
for i=1:1:ynum;
    for j=1:1:xnum;

         % Boundary conditions
        % Check distance from planetary centre
        dx=gridx(j)-xcentre;
        dy=gridy(i)-ycentre;
        dr=(dx^2+dy^2)^0.5;
        if (j==1 || j==xnum || i==1 || i==ynum || dr>=rbon)
            resfi(i,j)=0;
        else
        % Poisson equation d2FI/dx2+d2FI/dy2=KP*4*PI*GAMMA*RHO
            % Right Part
            resfi(i,j)=RP(i,j);
            % Left part
            % -dFI/dx(j-1/2)
            % Check distance from fi(i,j-1)to the planetary centre
             % Check distance from planetary centre
            dx1=gridx(j-1)-xcentre;
            dr1=(dx1^2+dy^2)^0.5;
            % fi(i,j-1) node is in the boundary condition domain
            if (dr1>=rbon)
                % Compute distance to the external boundary
                dx1=(rbon^2-dy^2)^0.5-abs(dx);
                % Add -dFI/dx derivative
                resfi(i,j)=resfi(i,j)+(finew(i,j)-vbon)/dx1/((xstp(j-1)+xstp(j))/2);
            else
            % fi(i,j-1) node is not in the boundary condition domain
                % Add -dFI/dx derivative
                resfi(i,j)=resfi(i,j)+(finew(i,j)-finew(i,j-1))/xstp(j-1)/((xstp(j-1)+xstp(j))/2); 
            end
            % dFI/dx(j+1/2)
            % Check distance from fi(i,j+1)to the planetary centre
             % Check distance from planetary centre
            dx1=gridx(j+1)-xcentre;
            dr1=(dx1^2+dy^2)^0.5;
            % fi(i,j+1) node is in the boundary condition domain
            if (dr1>=rbon)
                % Compute distance to the external boundary
                dx1=(rbon^2-dy^2)^0.5-abs(dx);
                % Add dFI/dx derivative 
                resfi(i,j)=resfi(i,j)-(vbon-finew(i,j))/dx1/((xstp(j-1)+xstp(j))/2);
            else
            % fi(i,j+1) node is not in the boundary condition domain
                % Add dFI/dx derivative to the left part
                resfi(i,j)=resfi(i,j)-(finew(i,j+1)-finew(i,j))/xstp(j)/((xstp(j-1)+xstp(j))/2); 
            end
            % -dFI/dy(i-1/2)
            % Check distance from fi(i-1,j)to the planetary centre
            % Check distance from planetary centre
            dy1=gridy(i-1)-ycentre;
            dr1=(dx^2+dy1^2)^0.5;
            % fi(i-1,j) node is in the boundary condition domain
            if (dr1>=rbon)
                % Compute distance to the external boundary
                dy1=(rbon^2-dx^2)^0.5-abs(dy);
                % Add dG/dx derivative to the left part
                resfi(i,j)=resfi(i,j)+(finew(i,j)-vbon)/dy1/((ystp(i-1)+ystp(i))/2);
            else
            % fi(i-1,j) node is not in the boundary condition domain
                % Add dFI/dx derivative to the left part
                resfi(i,j)=resfi(i,j)+(finew(i,j)-finew(i-1,j))/ystp(i-1)/((ystp(i-1)+ystp(i))/2);
            end
            % dFI/dy(i+1/2)
            % Check distance from fi(i+1,j)to the planetary centre
             % Check distance from planetary centre
            dy1=gridy(i+1)-ycentre;
            dr1=(dx^2+dy1^2)^0.5;
            % fi(i+1,j) node is in the boundary condition domain
            if (dr1>=rbon)
                % Compute distance to the external boundary
                dy1=(rbon^2-dx^2)^0.5-abs(dy);
                % Add dFI/dx derivative to the left part
                resfi(i,j)=resfi(i,j)-(vbon-finew(i,j))/dy1/((ystp(i-1)+ystp(i))/2);
            else
            % fi(i+1,j) node is not in the boundary condition domain
                % Add dFI/dx derivative to the left part
                resfi(i,j)=resfi(i,j)-(finew(i+1,j)-finew(i,j))/ystp(i)/((ystp(i-1)+ystp(i))/2);
            end
        end
    end
end


