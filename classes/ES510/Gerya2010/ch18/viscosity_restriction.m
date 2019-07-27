% Function Viscosity_restriction()
% This function makes restriction operation: 
% Interpolates shear (etas) and normal (etan) viscosity
% from finer (k) to coarcer (k+1) level
% using bilinear interpolation and given way of averaging (etamean)
% 0 - arithmetic
% 1 - geometric
% 2 - harmonic
% Resolution (xnum, ynum) and steps (xstp,ystp) for both levels
% is used for organizing the interpolation
function[etasc etanc]=Viscosity_restriction(etamean,k,xnum,ynum,xstp,ystp,etas,etan)
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
etasc=zeros(ynum(k+1),xnum(k+1));
etanc=zeros(ynum(k+1)-1,xnum(k+1)-1);
% Interpolation weigts
wts=zeros(ynum(k+1),xnum(k+1));
wtn=zeros(ynum(k+1)-1,xnum(k+1)-1);

% Interpolating viscosities from finer level nodes
% Cycle of node of finer (k) level
for i=1:1:ynum(k);
    for j=1:1:xnum(k);
        %  [ic,jc]---------+------[ic,jc+1]
        %     |            |          |
        %     |            |          |
        %     |----------(i,j)        |
        %     |                       |
        %  [ic+1,jc]----------------[ic+1,jc+1]
        
        % Shear viscosity (etas)
         % Defining horizontal and vertical positions of current vx(i,j) node
         % normalized to coarcer (k+1) grid steps
         xpos=(j-1)*xstp(k)/xstp(k+1);
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
         if (jc>xnum(k+1)-1)
             jc=xnum(k+1)-1;
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
         % Add shear viscosity from the finer level to the coarcer level 
         % for 4 nodes bounding the cell
         if (etamean==0)
            % Arithmetic mean
            etasc(ic,jc)=etasc(ic,jc)+(1.0-dx)*(1.0-dy)*etas(i,j);
            etasc(ic,jc+1)=etasc(ic,jc+1)+dx*(1.0-dy)*etas(i,j);
            etasc(ic+1,jc)=etasc(ic+1,jc)+(1.0-dx)*dy*etas(i,j);
            etasc(ic+1,jc+1)=etasc(ic+1,jc+1)+dx*dy*etas(i,j);
         else
            if (etamean==1)
                % Geometric mean
                etasc(ic,jc)=etasc(ic,jc)+(1.0-dx)*(1.0-dy)*log(etas(i,j));
                etasc(ic,jc+1)=etasc(ic,jc+1)+dx*(1.0-dy)*log(etas(i,j));
                etasc(ic+1,jc)=etasc(ic+1,jc)+(1.0-dx)*dy*log(etas(i,j));
                etasc(ic+1,jc+1)=etasc(ic+1,jc+1)+dx*dy*log(etas(i,j));
            else
                % Harmonic mean
                etasc(ic,jc)=etasc(ic,jc)+(1.0-dx)*(1.0-dy)/etas(i,j);
                etasc(ic,jc+1)=etasc(ic,jc+1)+dx*(1.0-dy)/etas(i,j);
                etasc(ic+1,jc)=etasc(ic+1,jc)+(1.0-dx)*dy/etas(i,j);
                etasc(ic+1,jc+1)=etasc(ic+1,jc+1)+dx*dy/etas(i,j);
            end
        end
        % Weights
        wts(ic,jc)=wts(ic,jc)+(1.0-dx)*(1.0-dy);
        wts(ic,jc+1)=wts(ic,jc+1)+dx*(1.0-dy);
        wts(ic+1,jc)=wts(ic+1,jc)+(1.0-dx)*dy;
        wts(ic+1,jc+1)=wts(ic+1,jc+1)+dx*dy;
        
        % Normal viscosity in P nodes
        if (i<ynum(k) && j<xnum(k));
            % Defining horizontal and vertical positions of current etan(i,j) node
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
            % Add normal viscosity from the finer level to the coarcer level 
            % for 4 nodes bounding the cell
            if (etamean==0)
                % Arithmetic mean
                etanc(ic,jc)=etanc(ic,jc)+(1.0-dx)*(1.0-dy)*etan(i,j);
                etanc(ic,jc+1)=etanc(ic,jc+1)+dx*(1.0-dy)*etan(i,j);
                etanc(ic+1,jc)=etanc(ic+1,jc)+(1.0-dx)*dy*etan(i,j);
                etanc(ic+1,jc+1)=etanc(ic+1,jc+1)+dx*dy*etan(i,j);
            else
                if (etamean==1)
                    % Geometric mean
                    etanc(ic,jc)=etanc(ic,jc)+(1.0-dx)*(1.0-dy)*log(etan(i,j));
                    etanc(ic,jc+1)=etanc(ic,jc+1)+dx*(1.0-dy)*log(etan(i,j));
                    etanc(ic+1,jc)=etanc(ic+1,jc)+(1.0-dx)*dy*log(etan(i,j));
                    etanc(ic+1,jc+1)=etanc(ic+1,jc+1)+dx*dy*log(etan(i,j));
                else
                    % Harmonic mean
                    etanc(ic,jc)=etanc(ic,jc)+(1.0-dx)*(1.0-dy)/etan(i,j);
                    etanc(ic,jc+1)=etanc(ic,jc+1)+dx*(1.0-dy)/etan(i,j);
                    etanc(ic+1,jc)=etanc(ic+1,jc)+(1.0-dx)*dy/etan(i,j);
                    etanc(ic+1,jc+1)=etanc(ic+1,jc+1)+dx*dy/etan(i,j);
                end
            end
            % Weights 
            wtn(ic,jc)=wtn(ic,jc)+(1.0-dx)*(1.0-dy);
            wtn(ic,jc+1)=wtn(ic,jc+1)+dx*(1.0-dy);
            wtn(ic+1,jc)=wtn(ic+1,jc)+(1.0-dx)*dy;
            wtn(ic+1,jc+1)=wtn(ic+1,jc+1)+dx*dy;
        end
        
    end            
end
% Recomputing viscosities
% for the coarcer level (k+1)
for ic=1:1:ynum(k+1);
    for jc=1:1:xnum(k+1);

        % Shear viscosity
        % Check for non-zero weights 
        if (wts(ic,jc)~=0);
            if (etamean==0)
                % Arithmetic mean
                etasc(ic,jc)=etasc(ic,jc)/wts(ic,jc);
            else
                if (etamean==1)
                    % Geometric mean
                    etasc(ic,jc)=exp(etasc(ic,jc)/wts(ic,jc));
                else
                    % Harmonic mean
                    etasc(ic,jc)=1/(etasc(ic,jc)/wts(ic,jc));
                end
            end
        end
        
        % Normal viscosity
        if (ic<ynum(k+1) && jc<xnum(k+1));
            if (wtn(ic,jc)~=0);
                if (etamean==0)
                    % Arithmetic mean
                    etanc(ic,jc)=etanc(ic,jc)/wtn(ic,jc);
                else
                    if (etamean==1)
                        % Geometric mean
                        etanc(ic,jc)=exp(etanc(ic,jc)/wtn(ic,jc));
                    else
                        % Harmonic mean
                        etanc(ic,jc)=1/(etanc(ic,jc)/wtn(ic,jc));
                    end
                end
            end
        end
        
    end            
end

