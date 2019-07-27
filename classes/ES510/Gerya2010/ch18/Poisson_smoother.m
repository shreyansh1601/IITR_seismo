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
function[phi,residual]=Poisson_smoother(iternum,krelax,xnum,ynum,xstp,ystp,R,phi)

% Poisson equation koefficients
xkf=1/xstp^2;
ykf=1/ystp^2;
xykf=-2/xstp^2-2/ystp^2;

% Gauss-Seidel iteration circle
for niter=1:1:iternum;
    % Solving of Poisson equation on nodes
    for i=1:1:ynum;
        for j=1:1:xnum;
            % Boundrary conditions 
            if (i==1 | i==ynum | j==1 | j==xnum);
                % Boundary condition for gravity potential
                phi(i,j)=0;
            % Solving Poisson equation
            else
                % Computing current residual
                residual(i,j)=R(i,j)-(xkf*(phi(i,j-1)+phi(i,j+1))+ykf*(phi(i-1,j)+phi(i+1,j))+xykf*phi(i,j));
                % Updating solution
                phi(i,j)=phi(i,j)+residual(i,j)/xykf*krelax;
            end
        end            
    end
    % End of Solving Poisson equation on nodes
end
% Computing final state of residuals
for i=1:1:ynum;
    for j=1:1:xnum;
        % Residuals for boundary conditions 
        if (i==1 | i==ynum | j==1 | j==xnum);
            residual(i,j)=0;
        % Solving Poisson equation
        else
        % Residuals for Poisson equation
            residual(i,j)=R(i,j)-(xkf*(phi(i,j-1)+phi(i,j+1))+ykf*(phi(i-1,j)+phi(i+1,j))+xykf*phi(i,j));
        end
    end            
end

