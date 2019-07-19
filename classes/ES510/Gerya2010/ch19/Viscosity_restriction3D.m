% Function Viscosity_restriction3D()
% This function makes restriction operation: 
% Interpolates shear (etaxy,etaxz,etayz) and normal (etan) viscosity
% from finer (n) to coarcer (n+1) level
% using bilinear interpolation and arithmetic averaging 
% 0 - arithmetic
% 1 - geometric
% 2 - harmonic
% Resolution (xnum, ynum) and steps (xstp,ystp) for both levels
% is used for organizing the interpolation
function[etaxyc etaxzc etayzc etanc]=Viscosity_restriction3D(n,xnum,ynum,znum,xstp,ystp,zstp,etaxy,etaxz,etayz,etan)
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
% Basic (density, shear viscosity - etas) nodes are shown with +
% normal viscosity (etan) is defined in pressure nodes P
% Ghost nodes shown outside the basic grid
% are used for boundary conditions

% Creating arrays for the coarser level
% Shear and normal viscosity
etaxyc=zeros(ynum(n+1),xnum(n+1),znum(n+1)-1);
etaxzc=zeros(ynum(n+1)-1,xnum(n+1),znum(n+1));
etayzc=zeros(ynum(n+1),xnum(n+1)-1,znum(n+1));
etanc=zeros(ynum(n+1)-1,xnum(n+1)-1,znum(n+1)-1);
% Interpolation weigts
wtxy=etaxyc;
wtxz=etaxzc;
wtyz=etayzc;
wtn=etanc;

% Interpolating viscosities from finer level nodes
% Cycle of node of finer (k) level
for i=1:1:ynum(n);
    for j=1:1:xnum(n);
        for k=1:1:znum(n);
            
            %  [ic,jc,kc]------+------[ic,jc+1,kc]
            %     |            |          |
            %     |            |          |
            %     |------------*          |
            %     |           /           |
            %     |       (i,j,k)         |
            %     |                       |
            %  [ic+1,jc,kc]-----------[ic+1,jc+1,kc]
        
        
            % Shear viscosity (etaxy)
            if (k<znum(n));
                % Defining horizontal and vertical positions of current etaxy(i,j,k) node
                % normalized to coarcer (n+1) grid steps
                xpos=(j-1)*xstp(n)/xstp(n+1);
                ypos=(i-1)*ystp(n)/ystp(n+1);
                zpos=(k-0.5)*zstp(n)/zstp(n+1)-0.5;
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
                if (jc>xnum(n+1)-1)
                    jc=xnum(n+1)-1;
                end
                if (ic<1)
                    ic=1;
                end
                if (ic>ynum(n+1)-1)
                    ic=ynum(n+1)-1;
                end
                if (kc<1)
                    kc=1;
                end
                if (kc>znum(n+1)-2)
                    kc=znum(n+1)-2;
                end
                % Define normalized distances from (i,j,k) node to [ic,jc,kc] node
                dx=xpos+1-jc;
                dy=ypos+1-ic;
                dz=zpos+1-kc;
                % Add shear viscosity (etaxy) from the finer level to the coarcer level 
                % for 8 nodes bounding the cell
                etaxyc(ic,jc,kc)=etaxyc(ic,jc,kc)+(1.0-dx)*(1.0-dy)*(1.0-dz)*etaxy(i,j,k);
                wtxy(ic,jc,kc)=wtxy(ic,jc,kc)+(1.0-dx)*(1.0-dy)*(1.0-dz);
                etaxyc(ic,jc,kc+1)=etaxyc(ic,jc,kc+1)+(1.0-dx)*(1.0-dy)*dz*etaxy(i,j,k);
                wtxy(ic,jc,kc+1)=wtxy(ic,jc,kc+1)+(1.0-dx)*(1.0-dy)*dz;
                etaxyc(ic+1,jc,kc)=etaxyc(ic+1,jc,kc)+(1.0-dx)*dy*(1.0-dz)*etaxy(i,j,k);
                wtxy(ic+1,jc,kc)=wtxy(ic+1,jc,kc)+(1.0-dx)*dy*(1.0-dz);
                etaxyc(ic+1,jc,kc+1)=etaxyc(ic+1,jc,kc+1)+(1.0-dx)*dy*dz*etaxy(i,j,k);
                wtxy(ic+1,jc,kc+1)=wtxy(ic+1,jc,kc+1)+(1.0-dx)*dy*dz;
                etaxyc(ic,jc+1,kc)=etaxyc(ic,jc+1,kc)+dx*(1.0-dy)*(1.0-dz)*etaxy(i,j,k);
                wtxy(ic,jc+1,kc)=wtxy(ic,jc+1,kc)+dx*(1.0-dy)*(1.0-dz);
                etaxyc(ic,jc+1,kc+1)=etaxyc(ic,jc+1,kc+1)+dx*(1.0-dy)*dz*etaxy(i,j,k);
                wtxy(ic,jc+1,kc+1)=wtxy(ic,jc+1,kc+1)+dx*(1.0-dy)*dz;
                etaxyc(ic+1,jc+1,kc)=etaxyc(ic+1,jc+1,kc)+dx*dy*(1.0-dz)*etaxy(i,j,k);
                wtxy(ic+1,jc+1,kc)=wtxy(ic+1,jc+1,kc)+dx*dy*(1.0-dz);
                etaxyc(ic+1,jc+1,kc+1)=etaxyc(ic+1,jc+1,kc+1)+dx*dy*dz*etaxy(i,j,k);
                wtxy(ic+1,jc+1,kc+1)=wtxy(ic+1,jc+1,kc+1)+dx*dy*dz;
            end
        
            % Shear viscosity (etaxz)
            if (i<ynum(n));
                % Defining horizontal and vertical positions of current etaxz(i,j,k) node
                % normalized to coarcer (n+1) grid steps
                xpos=(j-1)*xstp(n)/xstp(n+1);
                ypos=(i-0.5)*ystp(n)/ystp(n+1)-0.5;
                zpos=(k-1)*zstp(n)/zstp(n+1);
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
                if (jc>xnum(n+1)-1)
                    jc=xnum(n+1)-1;
                end
                if (ic<1)
                    ic=1;
                end
                if (ic>ynum(n+1)-2)
                    ic=ynum(n+1)-2;
                end
                if (kc<1)
                    kc=1;
                end
                if (kc>znum(n+1)-1)
                    kc=znum(n+1)-1;
                end
                % Define normalized distances from (i,j,k) node to [ic,jc,kc] node
                dx=xpos+1-jc;
                dy=ypos+1-ic;
                dz=zpos+1-kc;
                % Add shear viscosity (etaxy) from the finer level to the coarcer level 
                % for 8 nodes bounding the cell
                etaxzc(ic,jc,kc)=etaxzc(ic,jc,kc)+(1.0-dx)*(1.0-dy)*(1.0-dz)*etaxz(i,j,k);
                wtxz(ic,jc,kc)=wtxz(ic,jc,kc)+(1.0-dx)*(1.0-dy)*(1.0-dz);
                etaxzc(ic,jc,kc+1)=etaxzc(ic,jc,kc+1)+(1.0-dx)*(1.0-dy)*dz*etaxz(i,j,k);
                wtxz(ic,jc,kc+1)=wtxz(ic,jc,kc+1)+(1.0-dx)*(1.0-dy)*dz;
                etaxzc(ic+1,jc,kc)=etaxzc(ic+1,jc,kc)+(1.0-dx)*dy*(1.0-dz)*etaxz(i,j,k);
                wtxz(ic+1,jc,kc)=wtxz(ic+1,jc,kc)+(1.0-dx)*dy*(1.0-dz);
                etaxzc(ic+1,jc,kc+1)=etaxzc(ic+1,jc,kc+1)+(1.0-dx)*dy*dz*etaxz(i,j,k);
                wtxz(ic+1,jc,kc+1)=wtxz(ic+1,jc,kc+1)+(1.0-dx)*dy*dz;
                etaxzc(ic,jc+1,kc)=etaxzc(ic,jc+1,kc)+dx*(1.0-dy)*(1.0-dz)*etaxz(i,j,k);
                wtxz(ic,jc+1,kc)=wtxz(ic,jc+1,kc)+dx*(1.0-dy)*(1.0-dz);
                etaxzc(ic,jc+1,kc+1)=etaxzc(ic,jc+1,kc+1)+dx*(1.0-dy)*dz*etaxz(i,j,k);
                wtxz(ic,jc+1,kc+1)=wtxz(ic,jc+1,kc+1)+dx*(1.0-dy)*dz;
                etaxzc(ic+1,jc+1,kc)=etaxzc(ic+1,jc+1,kc)+dx*dy*(1.0-dz)*etaxz(i,j,k);
                wtxz(ic+1,jc+1,kc)=wtxz(ic+1,jc+1,kc)+dx*dy*(1.0-dz);
                etaxzc(ic+1,jc+1,kc+1)=etaxzc(ic+1,jc+1,kc+1)+dx*dy*dz*etaxz(i,j,k);
                wtxz(ic+1,jc+1,kc+1)=wtxz(ic+1,jc+1,kc+1)+dx*dy*dz;
            end

             % Shear viscosity (etayz)
            if (j<xnum(n));
                % Defining horizontal and vertical positions of current etaxz(i,j,k) node
                % normalized to coarcer (n+1) grid steps
                xpos=(j-0.5)*xstp(n)/xstp(n+1)-0.5;
                ypos=(i-1)*ystp(n)/ystp(n+1);
                zpos=(k-1)*zstp(n)/zstp(n+1);
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
                if (jc>xnum(n+1)-2)
                    jc=xnum(n+1)-2;
                end
                if (ic<1)
                    ic=1;
                end
                if (ic>ynum(n+1)-1)
                    ic=ynum(n+1)-1;
                end
                if (kc<1)
                    kc=1;
                end
                if (kc>znum(n+1)-1)
                    kc=znum(n+1)-1;
                end
                % Define normalized distances from (i,j,k) node to [ic,jc,kc] node
                dx=xpos+1-jc;
                dy=ypos+1-ic;
                dz=zpos+1-kc;
                % Add shear viscosity (etaxy) from the finer level to the coarcer level 
                % for 8 nodes bounding the cell
                etayzc(ic,jc,kc)=etayzc(ic,jc,kc)+(1.0-dx)*(1.0-dy)*(1.0-dz)*etayz(i,j,k);
                wtyz(ic,jc,kc)=wtyz(ic,jc,kc)+(1.0-dx)*(1.0-dy)*(1.0-dz);
                etayzc(ic,jc,kc+1)=etayzc(ic,jc,kc+1)+(1.0-dx)*(1.0-dy)*dz*etayz(i,j,k);
                wtyz(ic,jc,kc+1)=wtyz(ic,jc,kc+1)+(1.0-dx)*(1.0-dy)*dz;
                etayzc(ic+1,jc,kc)=etayzc(ic+1,jc,kc)+(1.0-dx)*dy*(1.0-dz)*etayz(i,j,k);
                wtyz(ic+1,jc,kc)=wtyz(ic+1,jc,kc)+(1.0-dx)*dy*(1.0-dz);
                etayzc(ic+1,jc,kc+1)=etayzc(ic+1,jc,kc+1)+(1.0-dx)*dy*dz*etayz(i,j,k);
                wtyz(ic+1,jc,kc+1)=wtyz(ic+1,jc,kc+1)+(1.0-dx)*dy*dz;
                etayzc(ic,jc+1,kc)=etayzc(ic,jc+1,kc)+dx*(1.0-dy)*(1.0-dz)*etayz(i,j,k);
                wtyz(ic,jc+1,kc)=wtyz(ic,jc+1,kc)+dx*(1.0-dy)*(1.0-dz);
                etayzc(ic,jc+1,kc+1)=etayzc(ic,jc+1,kc+1)+dx*(1.0-dy)*dz*etayz(i,j,k);
                wtyz(ic,jc+1,kc+1)=wtyz(ic,jc+1,kc+1)+dx*(1.0-dy)*dz;
                etayzc(ic+1,jc+1,kc)=etayzc(ic+1,jc+1,kc)+dx*dy*(1.0-dz)*etayz(i,j,k);
                wtyz(ic+1,jc+1,kc)=wtyz(ic+1,jc+1,kc)+dx*dy*(1.0-dz);
                etayzc(ic+1,jc+1,kc+1)=etayzc(ic+1,jc+1,kc+1)+dx*dy*dz*etayz(i,j,k);
                wtyz(ic+1,jc+1,kc+1)=wtyz(ic+1,jc+1,kc+1)+dx*dy*dz;
            end
           
            
            % Normal viscosity in P nodes
            if (i<ynum(n) && j<xnum(n) && k<znum(n));
                % Defining horizontal and vertical positions of current etan(i,j,k) node
                % normalized to coarcer (k+1) grid steps
                % Take intoaccount displacements in two staggered grids
                xpos=(j-0.5)*xstp(n)/xstp(n+1)-0.5;
                ypos=(i-0.5)*ystp(n)/ystp(n+1)-0.5;
                zpos=(k-0.5)*zstp(n)/zstp(n+1)-0.5;
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
                if (jc>xnum(n+1)-2)
                    jc=xnum(n+1)-2;
                end
                if (ic<1)
                    ic=1;
                end
                if (ic>ynum(n+1)-2)
                    ic=ynum(n+1)-2;
                end
                if (kc<1)
                    kc=1;
                end
                if (kc>znum(n+1)-2)
                    kc=znum(n+1)-2;
                end
                % Define normalized distances from (i,j) node to [ic,jc] node
                dx=xpos+1-jc;
                dy=ypos+1-ic;
                dz=zpos+1-kc;
                % Add normal viscosity from the finer level to the coarcer level 
                % for 8 nodes bounding the cell
                etanc(ic,jc,kc)=etanc(ic,jc,kc)+(1.0-dx)*(1.0-dy)*(1.0-dz)*etan(i,j,k);
                wtn(ic,jc,kc)=wtn(ic,jc,kc)+(1.0-dx)*(1.0-dy)*(1.0-dz);
                etanc(ic,jc,kc+1)=etanc(ic,jc,kc+1)+(1.0-dx)*(1.0-dy)*dz*etan(i,j,k);
                wtn(ic,jc,kc+1)=wtn(ic,jc,kc+1)+(1.0-dx)*(1.0-dy)*dz;
                etanc(ic+1,jc,kc)=etanc(ic+1,jc,kc)+(1.0-dx)*dy*(1.0-dz)*etan(i,j,k);
                wtn(ic+1,jc,kc)=wtn(ic+1,jc,kc)+(1.0-dx)*dy*(1.0-dz);
                etanc(ic+1,jc,kc+1)=etanc(ic+1,jc,kc+1)+(1.0-dx)*dy*dz*etan(i,j,k);
                wtn(ic+1,jc,kc+1)=wtn(ic+1,jc,kc+1)+(1.0-dx)*dy*dz;
                etanc(ic,jc+1,kc)=etanc(ic,jc+1,kc)+dx*(1.0-dy)*(1.0-dz)*etan(i,j,k);
                wtn(ic,jc+1,kc)=wtn(ic,jc+1,kc)+dx*(1.0-dy)*(1.0-dz);
                etanc(ic,jc+1,kc+1)=etanc(ic,jc+1,kc+1)+dx*(1.0-dy)*dz*etan(i,j,k);
                wtn(ic,jc+1,kc+1)=wtn(ic,jc+1,kc+1)+dx*(1.0-dy)*dz;
                etanc(ic+1,jc+1,kc)=etanc(ic+1,jc+1,kc)+dx*dy*(1.0-dz)*etan(i,j,k);
                wtn(ic+1,jc+1,kc)=wtn(ic+1,jc+1,kc)+dx*dy*(1.0-dz);
                etanc(ic+1,jc+1,kc+1)=etanc(ic+1,jc+1,kc+1)+dx*dy*dz*etan(i,j,k);
                wtn(ic+1,jc+1,kc+1)=wtn(ic+1,jc+1,kc+1)+dx*dy*dz;
            end

        end
    end            
end


% Recomputing viscosities
% for the coarcer level (k+1)
for ic=1:1:ynum(n+1);
    for jc=1:1:xnum(n+1);
        for kc=1:1:znum(n+1);

            % Shear viscosity (etaxy)
            if (kc<znum(n+1));
                % Check for non-zero weights 
                if (wtxy(ic,jc,kc)~=0);
                    etaxyc(ic,jc,kc)=etaxyc(ic,jc,kc)/wtxy(ic,jc,kc);
                end
            end
        
            % Shear viscosity (etaxz)
            if (ic<ynum(n+1));
                % Check for non-zero weights 
                if (wtxz(ic,jc,kc)~=0);
                    etaxzc(ic,jc,kc)=etaxzc(ic,jc,kc)/wtxz(ic,jc,kc);
                end
            end

            % Shear viscosity (etayz)
            if (jc<xnum(n+1));
                % Check for non-zero weights 
                if (wtyz(ic,jc,kc)~=0);
                    etayzc(ic,jc,kc)=etayzc(ic,jc,kc)/wtyz(ic,jc,kc);
                end
            end

            % Normal viscosity (etan)
            if (ic<ynum(n+1) && jc<xnum(n+1) && kc<znum(n+1));
                % Check for non-zero weights 
                if (wtn(ic,jc,kc)~=0);
                    etanc(ic,jc,kc)=etanc(ic,jc,kc)/wtn(ic,jc,kc);
                end
            end
            
        end
    end            
end

