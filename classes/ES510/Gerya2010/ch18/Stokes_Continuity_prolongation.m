% Function Stokes_prolongation()
% This function makes prolongation operation: 
% Interpolates corrections for solution (vx,vy, pr) from coarser (k) to finer (k-1) level
% using bilinear interpolation
% and produces updated solutions (dvx,dvy, dpr) for this finer level
%
% Resolution (xnum, ynum) and steps (xstp,ystp) for both levels
% is used for organizing the interpolation
function[dvx,dvy,dpr]=Stokes_prolongation(k,xnum,ynum,xstp,ystp,vx,vy,pr)
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

% Creating arrays for the finer level
% Solution updates
dvx=zeros(ynum(k-1)+1,xnum(k-1));
dvy=zeros(ynum(k-1),xnum(k-1)+1);
dpr=zeros(ynum(k-1)-1,xnum(k-1)-1);

% Interpolating residuals from finer level nodes
% Cycle of node of finer (k) level
for i=1:1:ynum(k-1)+1;
    for j=1:1:xnum(k-1)+1;
        %  [ic,jc]---------+------[ic,jc+1]
        %     |            |          |
        %     |            |          |
        %     |----------(i,j)        |
        %     |                       |
        %  [ic+1,jc]----------------[ic+1,jc+1]
        
        % x-Stokes equation residual
        if (j<xnum(k-1)+1);
            % Defining horizontal and vertical positions of current vx(i,j) node
            % normalized to coarcer (k+1) grid steps
            % Take intoaccount displacements in two staggered grids
            xpos=(j-1)*xstp(k-1)/xstp(k);
            ypos=(i-1.5)*ystp(k-1)/ystp(k)+0.5;
            % Defining smallest indexes [ic,jc] for the node 
            % bounding the cell on the coarcer grid 
            % in which current node of the finer grid is located
            % !!! SUBTRACT 0.5 since int16(0.5)=1
            jc=double(int16(xpos-0.5))+1;
            ic=double(int16(ypos-0.5))+1;
            % Check indexes
            if (jc<1)
                jc=1;
            end
            if (jc>xnum(k)-1)
                jc=xnum(k)-1;
            end
            if (ic<1)
              ic=1;
            end
            if (ic>ynum(k))
              ic=ynum(k);
            end
            % Define normalized distances from (i,j) node to [ic,jc] node
            dx=xpos+1-jc;
            dy=ypos+1-ic;
            % Add vx-velocity from the coarser to finer level 
            % for 4 nodes bounding the cell
            dvx(i,j)=dvx(i,j)+(1.0-dx)*(1.0-dy)*vx(ic,jc);
            dvx(i,j)=dvx(i,j)+dx*(1.0-dy)*vx(ic,jc+1);
            dvx(i,j)=dvx(i,j)+(1.0-dx)*dy*vx(ic+1,jc);
            dvx(i,j)=dvx(i,j)+dx*dy*vx(ic+1,jc+1);
        end
        
        % y-Stokes equation residual
        if (i<ynum(k-1)+1);
            % Defining horizontal and vertical positions of current vx(i,j) node
            % normalized to coarcer (k+1) grid steps
            % Take intoaccount displacements in two staggered grids
            xpos=(j-1.5)*xstp(k-1)/xstp(k)+0.5;
            ypos=(i-1)*ystp(k-1)/ystp(k);
            % Defining smallest indexes [ic,jc] for the node 
            % bounding the cell on the coarcer grid 
            % in which current node of the finer grid is located
            % !!! SUBTRACT 0.5 since int16(0.5)=1
            jc=double(int16(xpos-0.5))+1;
            ic=double(int16(ypos-0.5))+1;
            % Check indexes
            if (jc<1)
                jc=1;
            end
            if (jc>xnum(k))
                jc=xnum(k);
            end
            if (ic<1)
              ic=1;
            end
            if (ic>ynum(k)-1)
              ic=ynum(k)-1;
            end
            % Define normalized distances from (i,j) node to [ic,jc] node
            dx=xpos+1-jc;
            dy=ypos+1-ic;
            % Add vy-velocity from the coarser to finer level 
            % for 4 nodes bounding the cell
            dvy(i,j)=dvy(i,j)+(1.0-dx)*(1.0-dy)*vy(ic,jc);
            dvy(i,j)=dvy(i,j)+dx*(1.0-dy)*vy(ic,jc+1);
            dvy(i,j)=dvy(i,j)+(1.0-dx)*dy*vy(ic+1,jc);
            dvy(i,j)=dvy(i,j)+dx*dy*vy(ic+1,jc+1);
        end
        
        % Continuity equation residual
        if (i<ynum(k-1) && j<xnum(k-1));
            % Defining horizontal and vertical positions of current vx(i,j) node
            % normalized to coarcer (k+1) grid steps
            % Take intoaccount displacements in two staggered grids
            xpos=(j-0.5)*xstp(k-1)/xstp(k)-0.5;
            ypos=(i-0.5)*ystp(k-1)/ystp(k)-0.5;
            % Defining smallest indexes [ic,jc] for the node 
            % bounding the cell on the coarcer grid 
            % in which current node of the finer grid is located
            % !!! SUBTRACT 0.5 since int16(0.5)=1
            jc=double(int16(xpos-0.5))+1;
            ic=double(int16(ypos-0.5))+1;
            % Check indexes
            if (jc<1)
                jc=1;
            end
            if (jc>xnum(k)-2)
                jc=xnum(k)-2;
            end
            if (ic<1)
                ic=1;
            end
            if (ic>ynum(k)-2)
                ic=ynum(k)-2;
            end
            % Define normalized distances from (i,j) node to [ic,jc] node
            dx=xpos+1-jc;
            dy=ypos+1-ic;
            % Add pressure from the coarser to finer level 
            % for 4 nodes bounding the cell
            dpr(i,j)=dpr(i,j)+(1.0-dx)*(1.0-dy)*pr(ic,jc);
            dpr(i,j)=dpr(i,j)+dx*(1.0-dy)*pr(ic,jc+1);
            dpr(i,j)=dpr(i,j)+(1.0-dx)*dy*pr(ic+1,jc);
            dpr(i,j)=dpr(i,j)+dx*dy*pr(ic+1,jc+1);
        end
        
    end            
end
