% Function Stokes_Continuity_restriction()
% This function makes restriction operation: 
% Interpolates residuals (resx,resy,resc) from finer (k) to coarcer (k+1) level
% using bilinear interpolation
% and produces right parts (RX,RY,RC) for this coarcer level
%
% Resolution (xnum, ynum) and steps (xstp,ystp) for both levels
% is used for organizing the interpolation
function[RX,RY,RC]=Stokes_Continuity_restriction(k,xnum,ynum,xstp,ystp,resx,resy,resc)
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
RX=zeros(ynum(k+1)+1,xnum(k+1));
RY=zeros(ynum(k+1),xnum(k+1)+1);
RC=zeros(ynum(k+1)-1,xnum(k+1)-1);
% Interpolation weigts
wtx=zeros(ynum(k+1)+1,xnum(k+1));
wty=zeros(ynum(k+1),xnum(k+1)+1);
wtc=zeros(ynum(k+1)-1,xnum(k+1)-1);

% Interpolating residuals from finer level nodes
% Cycle of node of finer (k) level
for i=2:1:ynum(k);
    for j=2:1:xnum(k);
        %  [ic,jc]---------+------[ic,jc+1]
        %     |            |          |
        %     |            |          |
        %     |----------(i,j)        |
        %     |                       |
        %  [ic+1,jc]----------------[ic+1,jc+1]
        
        % x-Stokes equation residual
        if (j<xnum(k));
            % Defining horizontal and vertical positions of current vx(i,j) node
            % normalized to coarcer (k+1) grid steps
            % Take intoaccount displacements in two staggered grids
            xpos=(j-1)*xstp(k)/xstp(k+1);
            ypos=(i-1.5)*ystp(k)/ystp(k+1)+0.5;
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
            if (jc>xnum(k+1)-1)
                jc=xnum(k+1)-1;
            end
            if (ic<1)
              ic=1;
            end
            if (ic>ynum(k+1))
              ic=ynum(k+1);
            end
            % Define normalized distances from (i,j) node to [ic,jc] node
            dx=xpos+1-jc;
            dy=ypos+1-ic;
            % Add residual from the finer level to right parts of the coarcer level 
            % for 4 nodes bounding the cell
            RX(ic,jc)=RX(ic,jc)+(1.0-dx)*(1.0-dy)*resx(i,j);
            wtx(ic,jc)=wtx(ic,jc)+(1.0-dx)*(1.0-dy);
            RX(ic,jc+1)=RX(ic,jc+1)+dx*(1.0-dy)*resx(i,j);
            wtx(ic,jc+1)=wtx(ic,jc+1)+dx*(1.0-dy);
            RX(ic+1,jc)=RX(ic+1,jc)+(1.0-dx)*dy*resx(i,j);
            wtx(ic+1,jc)=wtx(ic+1,jc)+(1.0-dx)*dy;
            RX(ic+1,jc+1)=RX(ic+1,jc+1)+dx*dy*resx(i,j);
            wtx(ic+1,jc+1)=wtx(ic+1,jc+1)+dx*dy;
        end
        
        % y-Stokes equation residual
        if (i<ynum(k));
            % Defining horizontal and vertical positions of current vx(i,j) node
            % normalized to coarcer (k+1) grid steps
            % Take intoaccount displacements in two staggered grids
            xpos=(j-1.5)*xstp(k)/xstp(k+1)+0.5;
            ypos=(i-1)*ystp(k)/ystp(k+1);
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
            if (jc>xnum(k+1))
                jc=xnum(k+1);
            end
            if (ic<1)
              ic=1;
            end
            if (ic>ynum(k+1)-1)
              ic=ynum(k+1)-1;
            end
            % Define normalized distances from (i,j) node to [ic,jc] node
            dx=xpos+1-jc;
            dy=ypos+1-ic;
            % Add residual from the finer level to right parts of the coarcer level 
            % for 4 nodes bounding the cell
            RY(ic,jc)=RY(ic,jc)+(1.0-dx)*(1.0-dy)*resy(i,j);
            wty(ic,jc)=wty(ic,jc)+(1.0-dx)*(1.0-dy);
            RY(ic,jc+1)=RY(ic,jc+1)+dx*(1.0-dy)*resy(i,j);
            wty(ic,jc+1)=wty(ic,jc+1)+dx*(1.0-dy);
            RY(ic+1,jc)=RY(ic+1,jc)+(1.0-dx)*dy*resy(i,j);
            wty(ic+1,jc)=wty(ic+1,jc)+(1.0-dx)*dy;
            RY(ic+1,jc+1)=RY(ic+1,jc+1)+dx*dy*resy(i,j);
            wty(ic+1,jc+1)=wty(ic+1,jc+1)+dx*dy;
        end
        
        % Continuity equation residual
        if (i<ynum(k) && j<xnum(k));
            % Defining horizontal and vertical positions of current pr(i,j) node
            % normalized to coarcer (k+1) grid steps
            % Take intoaccount displacements in two staggered grids
            xpos=(j-0.5)*xstp(k)/xstp(k+1)-0.5;
            ypos=(i-0.5)*ystp(k)/ystp(k+1)-0.5;
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
            if (jc>xnum(k+1)-2)
            jc=xnum(k+1)-2;
            end
            if (ic<1)
            ic=1;
            end
            if (ic>ynum(k+1)-2)
                ic=ynum(k+1)-2;
            end
            % Define normalized distances from (i,j) node to [ic,jc] node
            dx=xpos+1-jc;
            dy=ypos+1-ic;
            % Add residual from the finer level to right parts of the coarcer level 
            % for 4 nodes bounding the cell
            RC(ic,jc)=RC(ic,jc)+(1.0-dx)*(1.0-dy)*resc(i,j);
            wtc(ic,jc)=wtc(ic,jc)+(1.0-dx)*(1.0-dy);
            RC(ic,jc+1)=RC(ic,jc+1)+dx*(1.0-dy)*resc(i,j);
            wtc(ic,jc+1)=wtc(ic,jc+1)+dx*(1.0-dy);
            RC(ic+1,jc)=RC(ic+1,jc)+(1.0-dx)*dy*resc(i,j);
            wtc(ic+1,jc)=wtc(ic+1,jc)+(1.0-dx)*dy;
            RC(ic+1,jc+1)=RC(ic+1,jc+1)+dx*dy*resc(i,j);
            wtc(ic+1,jc+1)=wtc(ic+1,jc+1)+dx*dy;
        end
        
    end            
end
% Recomputing right parts (RX, RY, RC)
% for the coarcer level (k+1)
for ic=1:1:ynum(k+1)+1;
    for jc=1:1:xnum(k+1)+1;

        % x-Stokes
        if (jc<xnum(k+1)+1);
            % Check for non-zero weights 
            if (wtx(ic,jc)~=0 && ic>1 && ic<ynum(k+1)+1 && jc>1 && jc<xnum(k+1));
                RX(ic,jc)=RX(ic,jc)/wtx(ic,jc);
            else
                RX(ic,jc)=0;
            end
        end
        
        % y-Stokes
        if (ic<ynum(k+1)+1);
            % Check for non-zero weights 
            if (wty(ic,jc)~=0 && ic>1 && ic<ynum(k+1) && jc>1 && jc<xnum(k+1)+1);
                RY(ic,jc)=RY(ic,jc)/wty(ic,jc);
            else
                RY(ic,jc)=0;
            end
        end

        % Continuity
        if (ic<ynum(k+1) && jc<xnum(k+1));
            if (wtc(ic,jc)~=0);
                RC(ic,jc)=RC(ic,jc)/wtc(ic,jc);
            else
                RC(ic,jc)=0;
            end
        end
        
    end            
end

