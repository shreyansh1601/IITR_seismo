% Function Stokes_Continuity3D_viscous_restriction()
% This function makes restriction operation: 
% Interpolates residuals (resx,resy,resc) from finer (n) to coarcer (n+1) level
% using bilinear interpolation
% and produces right parts (RX,RY,RZ,RC) for this coarcer level
% normal viscosity on finer (etanf) and coarser (etanc) levels
% is used for rescaling continuity residuals
%
% Resolution (xnum, ynum, znum) and steps (xstp,ystp,zstp) for both levels
% is used for organizing the interpolation
function[RX,RY,RZ,RC]=Stokes_Continuity3D_viscous_restriction(n,xnum,ynum,znum,xstp,ystp,zstp,resx,resy,resz,resc,etanf,etanc)
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

% Creating arrays for the coarser level
% Right parts
RX=zeros(ynum(n+1)+1,xnum(n+1),znum(n+1)+1);
RY=zeros(ynum(n+1),xnum(n+1)+1,znum(n+1)+1);
RZ=zeros(ynum(n+1)+1,xnum(n+1)+1,znum(n+1));
RC=zeros(ynum(n+1)-1,xnum(n+1)-1,znum(n+1)-1);
% Interpolation weigts
wtx=zeros(ynum(n+1)+1,xnum(n+1),znum(n+1)+1);
wty=zeros(ynum(n+1),xnum(n+1)+1,znum(n+1)+1);
wtz=zeros(ynum(n+1)+1,xnum(n+1)+1,znum(n+1));
wtc=zeros(ynum(n+1)-1,xnum(n+1)-1,znum(n+1)-1);

% Interpolating residuals from finer level nodes
% Cycle of node of finer (n) level
for i=2:1:ynum(n);
    for j=2:1:xnum(n);
        for k=2:1:znum(n);
          
            
            %  [ic,jc,kc]------+------[ic,jc+1,kc]
            %     |            |          |
            %     |            |          |
            %     |------------*          |
            %     |           /           |
            %     |       (i,j,k)         |
            %     |                       |
            %  [ic+1,jc,kc]-----------[ic+1,jc+1,kc]
        
            % x-Stokes equation residual
            if (j<xnum(n));
                % Defining horizontal and vertical positions of current vx(i,j,k) node
                % normalized to coarcer (n+1) grid steps
                % Take into account displacements in two staggered grids
                xpos=(j-1)*xstp(n)/xstp(n+1);
                ypos=(i-1.5)*ystp(n)/ystp(n+1)+0.5;
                zpos=(k-1.5)*zstp(n)/zstp(n+1)+0.5;
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
                if (ic>ynum(n+1))
                    ic=ynum(n+1);
                end
                if (kc<1)
                    kc=1;
                end
                if (kc>znum(n+1))
                    kc=znum(n+1);
                end
                % Define normalized distances from (i,j,k) node to [ic,jc,kc] node
                dx=xpos+1-jc;
                dy=ypos+1-ic;
                dz=zpos+1-kc;
                % Add residual from the finer level to right parts of the coarcer level 
                % for 8 nodes bounding the cell
                RX(ic,jc,kc)=RX(ic,jc,kc)+(1.0-dx)*(1.0-dy)*(1.0-dz)*resx(i,j,k);
                wtx(ic,jc,kc)=wtx(ic,jc,kc)+(1.0-dx)*(1.0-dy)*(1.0-dz);
                RX(ic,jc,kc+1)=RX(ic,jc,kc+1)+(1.0-dx)*(1.0-dy)*dz*resx(i,j,k);
                wtx(ic,jc,kc+1)=wtx(ic,jc,kc+1)+(1.0-dx)*(1.0-dy)*dz;
                RX(ic+1,jc,kc)=RX(ic+1,jc,kc)+(1.0-dx)*dy*(1.0-dz)*resx(i,j,k);
                wtx(ic+1,jc,kc)=wtx(ic+1,jc,kc)+(1.0-dx)*dy*(1.0-dz);
                RX(ic+1,jc,kc+1)=RX(ic+1,jc,kc+1)+(1.0-dx)*dy*dz*resx(i,j,k);
                wtx(ic+1,jc,kc+1)=wtx(ic+1,jc,kc+1)+(1.0-dx)*dy*dz;
                RX(ic,jc+1,kc)=RX(ic,jc+1,kc)+dx*(1.0-dy)*(1.0-dz)*resx(i,j,k);
                wtx(ic,jc+1,kc)=wtx(ic,jc+1,kc)+dx*(1.0-dy)*(1.0-dz);
                RX(ic,jc+1,kc+1)=RX(ic,jc+1,kc+1)+dx*(1.0-dy)*dz*resx(i,j,k);
                wtx(ic,jc+1,kc+1)=wtx(ic,jc+1,kc+1)+dx*(1.0-dy)*dz;
                RX(ic+1,jc+1,kc)=RX(ic+1,jc+1,kc)+dx*dy*(1.0-dz)*resx(i,j,k);
                wtx(ic+1,jc+1,kc)=wtx(ic+1,jc+1,kc)+dx*dy*(1.0-dz);
                RX(ic+1,jc+1,kc+1)=RX(ic+1,jc+1,kc+1)+dx*dy*dz*resx(i,j,k);
                wtx(ic+1,jc+1,kc+1)=wtx(ic+1,jc+1,kc+1)+dx*dy*dz;
            end
        
            % y-Stokes equation residual
            if (i<ynum(n));
                % Defining horizontal and vertical positions of current vy(i,j,k) node
                % normalized to coarcer (n+1) grid steps
                % Take into account displacements in two staggered grids
                xpos=(j-1.5)*xstp(n)/xstp(n+1)+0.5;
                ypos=(i-1)*ystp(n)/ystp(n+1);
                zpos=(k-1.5)*zstp(n)/zstp(n+1)+0.5;
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
                if (jc>xnum(n+1))
                    jc=xnum(n+1);
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
                if (kc>znum(n+1))
                    kc=znum(n+1);
                end
                % Define normalized distances from (i,j,k) node to [ic,jc,kc] node
                dx=xpos+1-jc;
                dy=ypos+1-ic;
                dz=zpos+1-kc;
                % Add residual from the finer level to right parts of the coarcer level 
                % for 8 nodes bounding the cell
                RY(ic,jc,kc)=RY(ic,jc,kc)+(1.0-dx)*(1.0-dy)*(1.0-dz)*resy(i,j,k);
                wty(ic,jc,kc)=wty(ic,jc,kc)+(1.0-dx)*(1.0-dy)*(1.0-dz);
                RY(ic,jc,kc+1)=RY(ic,jc,kc+1)+(1.0-dx)*(1.0-dy)*dz*resy(i,j,k);
                wty(ic,jc,kc+1)=wty(ic,jc,kc+1)+(1.0-dx)*(1.0-dy)*dz;
                RY(ic+1,jc,kc)=RY(ic+1,jc,kc)+(1.0-dx)*dy*(1.0-dz)*resy(i,j,k);
                wty(ic+1,jc,kc)=wty(ic+1,jc,kc)+(1.0-dx)*dy*(1.0-dz);
                RY(ic+1,jc,kc+1)=RY(ic+1,jc,kc+1)+(1.0-dx)*dy*dz*resy(i,j,k);
                wty(ic+1,jc,kc+1)=wty(ic+1,jc,kc+1)+(1.0-dx)*dy*dz;
                RY(ic,jc+1,kc)=RY(ic,jc+1,kc)+dx*(1.0-dy)*(1.0-dz)*resy(i,j,k);
                wty(ic,jc+1,kc)=wty(ic,jc+1,kc)+dx*(1.0-dy)*(1.0-dz);
                RY(ic,jc+1,kc+1)=RY(ic,jc+1,kc+1)+dx*(1.0-dy)*dz*resy(i,j,k);
                wty(ic,jc+1,kc+1)=wty(ic,jc+1,kc+1)+dx*(1.0-dy)*dz;
                RY(ic+1,jc+1,kc)=RY(ic+1,jc+1,kc)+dx*dy*(1.0-dz)*resy(i,j,k);
                wty(ic+1,jc+1,kc)=wty(ic+1,jc+1,kc)+dx*dy*(1.0-dz);
                RY(ic+1,jc+1,kc+1)=RY(ic+1,jc+1,kc+1)+dx*dy*dz*resy(i,j,k);
                wty(ic+1,jc+1,kc+1)=wty(ic+1,jc+1,kc+1)+dx*dy*dz;
            end
        
            % z-Stokes equation residual
            if (k<znum(n));
                % Defining horizontal and vertical positions of current vz(i,j,k) node
                % normalized to coarcer (n+1) grid steps
                % Take into account displacements in two staggered grids
                xpos=(j-1.5)*xstp(n)/xstp(n+1)+0.5;
                ypos=(i-1.5)*ystp(n)/ystp(n+1)+0.5;
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
                if (jc>xnum(n+1))
                    jc=xnum(n+1);
                end
                if (ic<1)
                    ic=1;
                end
                if (ic>ynum(n+1))
                    ic=ynum(n+1);
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
                % Add residual from the finer level to right parts of the coarcer level 
                % for 8 nodes bounding the cell
                RZ(ic,jc,kc)=RZ(ic,jc,kc)+(1.0-dx)*(1.0-dy)*(1.0-dz)*resz(i,j,k);
                wtz(ic,jc,kc)=wtz(ic,jc,kc)+(1.0-dx)*(1.0-dy)*(1.0-dz);
                RZ(ic,jc,kc+1)=RZ(ic,jc,kc+1)+(1.0-dx)*(1.0-dy)*dz*resz(i,j,k);
                wtz(ic,jc,kc+1)=wtz(ic,jc,kc+1)+(1.0-dx)*(1.0-dy)*dz;
                RZ(ic+1,jc,kc)=RZ(ic+1,jc,kc)+(1.0-dx)*dy*(1.0-dz)*resz(i,j,k);
                wtz(ic+1,jc,kc)=wtz(ic+1,jc,kc)+(1.0-dx)*dy*(1.0-dz);
                RZ(ic+1,jc,kc+1)=RZ(ic+1,jc,kc+1)+(1.0-dx)*dy*dz*resz(i,j,k);
                wtz(ic+1,jc,kc+1)=wtz(ic+1,jc,kc+1)+(1.0-dx)*dy*dz;
                RZ(ic,jc+1,kc)=RZ(ic,jc+1,kc)+dx*(1.0-dy)*(1.0-dz)*resz(i,j,k);
                wtz(ic,jc+1,kc)=wtz(ic,jc+1,kc)+dx*(1.0-dy)*(1.0-dz);
                RZ(ic,jc+1,kc+1)=RZ(ic,jc+1,kc+1)+dx*(1.0-dy)*dz*resz(i,j,k);
                wtz(ic,jc+1,kc+1)=wtz(ic,jc+1,kc+1)+dx*(1.0-dy)*dz;
                RZ(ic+1,jc+1,kc)=RZ(ic+1,jc+1,kc)+dx*dy*(1.0-dz)*resz(i,j,k);
                wtz(ic+1,jc+1,kc)=wtz(ic+1,jc+1,kc)+dx*dy*(1.0-dz);
                RZ(ic+1,jc+1,kc+1)=RZ(ic+1,jc+1,kc+1)+dx*dy*dz*resz(i,j,k);
                wtz(ic+1,jc+1,kc+1)=wtz(ic+1,jc+1,kc+1)+dx*dy*dz;
            end
        
            % Continuity equation residual
            if (i<ynum(n) && j<xnum(n) && k<znum(n));
                % Defining horizontal and vertical positions of current pr(i,j,k) node
                % normalized to coarcer (n+1) grid steps
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
                % Define normalized distances from (i,j,k) node to [ic,jc] node
                dx=xpos+1-jc;
                dy=ypos+1-ic;
                dz=zpos+1-kc;
                % Add residual from the finer level to right parts of the coarcer level 
                % for 8 nodes bounding the cell
                RC(ic,jc,kc)=RC(ic,jc,kc)+(1.0-dx)*(1.0-dy)*(1.0-dz)*resc(i,j,k)*etanf(i,j,k);
                wtc(ic,jc,kc)=wtc(ic,jc,kc)+(1.0-dx)*(1.0-dy)*(1.0-dz);
                RC(ic,jc,kc+1)=RC(ic,jc,kc+1)+(1.0-dx)*(1.0-dy)*dz*resc(i,j,k)*etanf(i,j,k);
                wtc(ic,jc,kc+1)=wtc(ic,jc,kc+1)+(1.0-dx)*(1.0-dy)*dz;
                RC(ic+1,jc,kc)=RC(ic+1,jc,kc)+(1.0-dx)*dy*(1.0-dz)*resc(i,j,k)*etanf(i,j,k);
                wtc(ic+1,jc,kc)=wtc(ic+1,jc,kc)+(1.0-dx)*dy*(1.0-dz);
                RC(ic+1,jc,kc+1)=RC(ic+1,jc,kc+1)+(1.0-dx)*dy*dz*resc(i,j,k)*etanf(i,j,k);
                wtc(ic+1,jc,kc+1)=wtc(ic+1,jc,kc+1)+(1.0-dx)*dy*dz;
                RC(ic,jc+1,kc)=RC(ic,jc+1,kc)+dx*(1.0-dy)*(1.0-dz)*resc(i,j,k)*etanf(i,j,k);
                wtc(ic,jc+1,kc)=wtc(ic,jc+1,kc)+dx*(1.0-dy)*(1.0-dz);
                RC(ic,jc+1,kc+1)=RC(ic,jc+1,kc+1)+dx*(1.0-dy)*dz*resc(i,j,k)*etanf(i,j,k);
                wtc(ic,jc+1,kc+1)=wtc(ic,jc+1,kc+1)+dx*(1.0-dy)*dz;
                RC(ic+1,jc+1,kc)=RC(ic+1,jc+1,kc)+dx*dy*(1.0-dz)*resc(i,j,k)*etanf(i,j,k);
                wtc(ic+1,jc+1,kc)=wtc(ic+1,jc+1,kc)+dx*dy*(1.0-dz);
                RC(ic+1,jc+1,kc+1)=RC(ic+1,jc+1,kc+1)+dx*dy*dz*resc(i,j,k)*etanf(i,j,k);
                wtc(ic+1,jc+1,kc+1)=wtc(ic+1,jc+1,kc+1)+dx*dy*dz;
            end
        
        end
    end            
end
% Recomputing right parts (RX, RY, RZ, RC)
% for the coarcer level (n+1)
for ic=1:1:ynum(n+1)+1;
    for jc=1:1:xnum(n+1)+1;
        for kc=1:1:znum(n+1)+1;

            % x-Stokes
            if (jc<xnum(n+1)+1);
                % Check for non-zero weights 
                if (wtx(ic,jc,kc)~=0 && ic>1 && ic<ynum(n+1)+1 && jc>1 && jc<xnum(n+1) && kc>1 && kc<znum(n+1)+1);
                    RX(ic,jc,kc)=RX(ic,jc,kc)/wtx(ic,jc,kc);
                else
                    RX(ic,jc,kc)=0;
                end
            end
        
            % y-Stokes
            if (ic<ynum(n+1)+1);
                % Check for non-zero weights 
                if (wty(ic,jc,kc)~=0 && ic>1 && ic<ynum(n+1) && jc>1 && jc<xnum(n+1)+1 && kc>1 && kc<znum(n+1)+1);
                    RY(ic,jc,kc)=RY(ic,jc,kc)/wty(ic,jc,kc);
                else
                    RY(ic,jc,kc)=0;
                end
            end

            % z-Stokes
            if (kc<znum(n+1)+1);
                % Check for non-zero weights 
                if (wtz(ic,jc,kc)~=0 && ic>1 && ic<ynum(n+1)+1 && jc>1 && jc<xnum(n+1)+1 && kc>1 && kc<znum(n+1));
                    RZ(ic,jc,kc)=RZ(ic,jc,kc)/wtz(ic,jc,kc);
                else
                    RZ(ic,jc,kc)=0;
                end
            end

            % Continuity
            if (ic<ynum(n+1) && jc<xnum(n+1) && kc<znum(n+1));
                if (wtc(ic,jc,kc)~=0);
                    RC(ic,jc,kc)=RC(ic,jc,kc)/wtc(ic,jc,kc)/etanc(ic,jc,kc);
                else
                    RC(ic,jc,kc)=0;
                end
            end
        
        end
    end            
end

