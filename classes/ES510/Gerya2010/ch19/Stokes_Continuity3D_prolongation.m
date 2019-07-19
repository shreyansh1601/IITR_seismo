% Function Stokes_Continuity3D_prolongation()
% This function makes prolongation operation: 
% Interpolates corrections for solution (vx,vy,vz,pr) from coarser (n) to finer (n-1) level
% using bilinear interpolation
% and produces updated solutions (dvx,dvy,dvz,dpr) for this finer level
%
% Resolution (xnum,ynum,znum) and steps (xstp,ystp,zstp) for both levels
% is used for organizing the interpolation
function[dvx,dvy,dvz,dpr]=Stokes_Continuity3D_prolongation(n,xnum,ynum,znum,xstp,ystp,zstp,vx,vy,vz,pr)
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

% Creating arrays for the finer (n-1) level
% Solution updates
dvx=zeros(ynum(n-1)+1,xnum(n-1),znum(n-1)+1);
dvy=zeros(ynum(n-1),xnum(n-1)+1,znum(n-1)+1);
dvz=zeros(ynum(n-1)+1,xnum(n-1)+1,znum(n-1));
dpr=zeros(ynum(n-1)-1,xnum(n-1)-1,znum(n-1)-1);

% Interpolating solutions from coarser level nodes
% Cycle of node of finer (k) level
for i=1:1:ynum(n-1)+1;
    for j=1:1:xnum(n-1)+1;
        for k=1:1:znum(n-1)+1;
         
            
            %  [ic,jc,kc]------+------[ic,jc+1,kc]
            %     |            |          |
            %     |            |          |
            %     |------------*          |
            %     |           /           |
            %     |       (i,j,k)         |
            %     |                       |
            %  [ic+1,jc,kc]-----------[ic+1,jc+1,kc]
        
            % x-Stokes equation residual
            if (j<xnum(n-1)+1);
                % Defining horizontal and vertical positions of current vx(i,j,k) node
                % normalized to coarcer (n) grid steps
                % Take intoaccount displacements in two staggered grids
                xpos=(j-1)*xstp(n-1)/xstp(n);
                ypos=(i-1.5)*ystp(n-1)/ystp(n)+0.5;
                zpos=(k-1.5)*zstp(n-1)/zstp(n)+0.5;
                % Defining smallest indexes [ic,jc,kc] for the node 
                % bounding the cell on the coarcer grid 
                % in which current node of the finer grid is located
                % !!! SUBTRACT 0.5 since int16(0.5)=1
                jc=double(int16(xpos-0.5))+1;
                ic=double(int16(ypos-0.5))+1;
                kc=double(int16(zpos-0.5))+1;
                % Check indexes
                if (jc<1)
                    jc=1;
                end
                if (jc>xnum(n)-1)
                    jc=xnum(n)-1;
                end
                if (ic<1)
                    ic=1;
                end
                if (ic>ynum(n))
                    ic=ynum(n);
                end
                if (kc<1)
                    kc=1;
                end
                if (kc>znum(n))
                    kc=znum(n);
                end
                % Define normalized distances from (i,j,k) node to [ic,jc,kc] node
                dx=xpos+1-jc;
                dy=ypos+1-ic;
                dz=zpos+1-kc;
                % Add vx-velocity from the coarser to finer level 
                % using bilinear interpolation from 8 nodes bounding the cell
                dvx(i,j,k)=dvx(i,j,k)+(1.0-dx)*(1.0-dy)*(1.0-dz)*vx(ic,jc,kc);
                dvx(i,j,k)=dvx(i,j,k)+dx*(1.0-dy)*(1.0-dz)*vx(ic,jc+1,kc);
                dvx(i,j,k)=dvx(i,j,k)+(1.0-dx)*(1.0-dz)*dy*vx(ic+1,jc,kc);
                dvx(i,j,k)=dvx(i,j,k)+dx*dy*(1.0-dz)*vx(ic+1,jc+1,kc);
                dvx(i,j,k)=dvx(i,j,k)+(1.0-dx)*(1.0-dy)*dz*vx(ic,jc,kc+1);
                dvx(i,j,k)=dvx(i,j,k)+dx*(1.0-dy)*dz*vx(ic,jc+1,kc+1);
                dvx(i,j,k)=dvx(i,j,k)+(1.0-dx)*dz*dy*vx(ic+1,jc,kc+1);
                dvx(i,j,k)=dvx(i,j,k)+dx*dy*dz*vx(ic+1,jc+1,kc+1);
            end
        
            % y-Stokes equation residual
            if (i<ynum(n-1)+1);
                % Defining horizontal and vertical positions of current vy(i,j,k) node
                % normalized to coarcer (n) grid steps
                % Take intoaccount displacements in two staggered grids
                xpos=(j-1.5)*xstp(n-1)/xstp(n)+0.5;
                ypos=(i-1)*ystp(n-1)/ystp(n);
                zpos=(k-1.5)*zstp(n-1)/zstp(n)+0.5;
                % Defining smallest indexes [ic,jc,kc] for the node 
                % bounding the cell on the coarcer grid 
                % in which current node of the finer grid is located
                % !!! SUBTRACT 0.5 since int16(0.5)=1
                jc=double(int16(xpos-0.5))+1;
                ic=double(int16(ypos-0.5))+1;
                kc=double(int16(zpos-0.5))+1;
                % Check indexes
                if (jc<1)
                    jc=1;
                end
                if (jc>xnum(n))
                    jc=xnum(n);
                end
                if (ic<1)
                    ic=1;
                end
                if (ic>ynum(n)-1)
                    ic=ynum(n)-1;
                end
                if (kc<1)
                    kc=1;
                end
                if (kc>znum(n))
                    kc=znum(n);
                end
                % Define normalized distances from (i,j,k) node to [ic,jc,kc] node
                dx=xpos+1-jc;
                dy=ypos+1-ic;
                dz=zpos+1-kc;
                % Add vy-velocity from the coarser to finer level 
                % using bilinear interpolation from 8 nodes bounding the cell
                dvy(i,j,k)=dvy(i,j,k)+(1.0-dx)*(1.0-dy)*(1.0-dz)*vy(ic,jc,kc);
                dvy(i,j,k)=dvy(i,j,k)+dx*(1.0-dy)*(1.0-dz)*vy(ic,jc+1,kc);
                dvy(i,j,k)=dvy(i,j,k)+(1.0-dx)*(1.0-dz)*dy*vy(ic+1,jc,kc);
                dvy(i,j,k)=dvy(i,j,k)+dx*dy*(1.0-dz)*vy(ic+1,jc+1,kc);
                dvy(i,j,k)=dvy(i,j,k)+(1.0-dx)*(1.0-dy)*dz*vy(ic,jc,kc+1);
                dvy(i,j,k)=dvy(i,j,k)+dx*(1.0-dy)*dz*vy(ic,jc+1,kc+1);
                dvy(i,j,k)=dvy(i,j,k)+(1.0-dx)*dz*dy*vy(ic+1,jc,kc+1);
                dvy(i,j,k)=dvy(i,j,k)+dx*dy*dz*vy(ic+1,jc+1,kc+1);
            end

            % z-Stokes equation residual
            if (k<znum(n-1)+1);
                % Defining horizontal and vertical positions of current vx(i,j,k) node
                % normalized to coarcer (n) grid steps
                % Take intoaccount displacements in two staggered grids
                xpos=(j-1.5)*xstp(n-1)/xstp(n)+0.5;
                ypos=(i-1.5)*ystp(n-1)/ystp(n)+0.5;
                zpos=(k-1)*zstp(n-1)/zstp(n);
                % Defining smallest indexes [ic,jc,kc] for the node 
                % bounding the cell on the coarcer grid 
                % in which current node of the finer grid is located
                % !!! SUBTRACT 0.5 since int16(0.5)=1
                jc=double(int16(xpos-0.5))+1;
                ic=double(int16(ypos-0.5))+1;
                kc=double(int16(zpos-0.5))+1;
                % Check indexes
                if (jc<1)
                    jc=1;
                end
                if (jc>xnum(n))
                    jc=xnum(n);
                end
                if (ic<1)
                    ic=1;
                end
                if (ic>ynum(n))
                    ic=ynum(n);
                end
                if (kc<1)
                    kc=1;
                end
                if (kc>znum(n)-1)
                    kc=znum(n)-1;
                end
                % Define normalized distances from (i,j,k) node to [ic,jc,kc] node
                dx=xpos+1-jc;
                dy=ypos+1-ic;
                dz=zpos+1-kc;
                % Add vz-velocity from the coarser to finer level 
                % using bilinear interpolation from 8 nodes bounding the cell
                dvz(i,j,k)=dvz(i,j,k)+(1.0-dx)*(1.0-dy)*(1.0-dz)*vz(ic,jc,kc);
                dvz(i,j,k)=dvz(i,j,k)+dx*(1.0-dy)*(1.0-dz)*vz(ic,jc+1,kc);
                dvz(i,j,k)=dvz(i,j,k)+(1.0-dx)*(1.0-dz)*dy*vz(ic+1,jc,kc);
                dvz(i,j,k)=dvz(i,j,k)+dx*dy*(1.0-dz)*vz(ic+1,jc+1,kc);
                dvz(i,j,k)=dvz(i,j,k)+(1.0-dx)*(1.0-dy)*dz*vz(ic,jc,kc+1);
                dvz(i,j,k)=dvz(i,j,k)+dx*(1.0-dy)*dz*vz(ic,jc+1,kc+1);
                dvz(i,j,k)=dvz(i,j,k)+(1.0-dx)*dz*dy*vz(ic+1,jc,kc+1);
                dvz(i,j,k)=dvz(i,j,k)+dx*dy*dz*vz(ic+1,jc+1,kc+1);
            end
        
            % Continuity equation residual
            if (i<ynum(n-1) && j<xnum(n-1) && k<znum(n-1));
                % Defining horizontal and vertical positions of current vx(i,j,k) node
                % normalized to coarcer (n+1) grid steps
                % Take intoaccount displacements in two staggered grids
                xpos=(j-0.5)*xstp(n-1)/xstp(n)-0.5;
                ypos=(i-0.5)*ystp(n-1)/ystp(n)-0.5;
                zpos=(k-0.5)*zstp(n-1)/zstp(n)-0.5;
                % Defining smallest indexes [ic,jc,kc] for the node 
                % bounding the cell on the coarcer grid 
                % in which current node of the finer grid is located
                % !!! SUBTRACT 0.5 since int16(0.5)=1
                jc=double(int16(xpos-0.5))+1;
                ic=double(int16(ypos-0.5))+1;
                kc=double(int16(zpos-0.5))+1;
                % Check indexes
                if (jc<1)
                    jc=1;
                end
                if (jc>xnum(n)-2)
                    jc=xnum(n)-2;
                end
                if (ic<1)
                    ic=1;
                end
                if (ic>ynum(n)-2)
                    ic=ynum(n)-2;
                end
                if (kc<1)
                    kc=1;
                end
                if (kc>znum(n)-2)
                    kc=znum(n)-2;
                end
                % Define normalized distances from (i,j,k) node to [ic,jc,kc] node
                dx=xpos+1-jc;
                dy=ypos+1-ic;
                dz=zpos+1-kc;
                % Add pressure from the coarser to finer level 
                % using bilinear interpolation from 8 nodes bounding the cell
                dpr(i,j,k)=dpr(i,j,k)+(1.0-dx)*(1.0-dy)*(1.0-dz)*pr(ic,jc,kc);
                dpr(i,j,k)=dpr(i,j,k)+dx*(1.0-dy)*(1.0-dz)*pr(ic,jc+1,kc);
                dpr(i,j,k)=dpr(i,j,k)+(1.0-dx)*(1.0-dz)*dy*pr(ic+1,jc,kc);
                dpr(i,j,k)=dpr(i,j,k)+dx*dy*(1.0-dz)*pr(ic+1,jc+1,kc);
                dpr(i,j,k)=dpr(i,j,k)+(1.0-dx)*(1.0-dy)*dz*pr(ic,jc,kc+1);
                dpr(i,j,k)=dpr(i,j,k)+dx*(1.0-dy)*dz*pr(ic,jc+1,kc+1);
                dpr(i,j,k)=dpr(i,j,k)+(1.0-dx)*dz*dy*pr(ic+1,jc,kc+1);
                dpr(i,j,k)=dpr(i,j,k)+dx*dy*dz*pr(ic+1,jc+1,kc+1);
            end
        
        end
    end            
end

