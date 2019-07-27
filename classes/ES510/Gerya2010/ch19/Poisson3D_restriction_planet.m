% Function Poisson3D_restriction_planet()
% This function makes restriction operation: 
% Interpolates residuals (residual) from finer (n) to coarcer (n+1) level
% using bilinear interpolation
% and produces right parts (R) for this coarcer level
%
% Resolution (xnum, ynum) and steps (xstp,ystp,zstp) for both levels
% is used for organizing the interpolation
function[R]=Poisson3D_restriction_planet(n,xnum,ynum,znum,xstp,ystp,zstp,residual,bonf,bonc)

% Creating arrays for the coarser level
% Right parts
R=zeros(ynum(n+1),xnum(n+1),znum(n+1));
% Interpolation weigts
wt=zeros(ynum(n+1),xnum(n+1),znum(n+1));

% Interpolating residuals from finer level nodes
% Cycle of node of finer (k) level
for i=2:1:ynum(n)-1
    for j=2:1:xnum(n)-1
        for k=2:1:znum(n)-1
            % Skip nodes outside of the boundary circle
            if(bonf(i,j,k)>0)
                % Defining horizontal and vertical positions of current (i,j) node
                % normalized to coarcer (k+1) grid steps
                xpos=(j-1)*xstp(n)/xstp(n+1);
                ypos=(i-1)*ystp(n)/ystp(n+1);
                zpos=(k-1)*zstp(n)/zstp(n+1);
        
                %  [ic,jc,kc]------+------[ic,jc+1,kc]
                %     |            |          |
                %     |            |          |
                %     |------------*          |
                %     |           /           |
                %     |       (i,j,k)         |
                %     |                       |
                %  [ic+1,jc,kc]-----------[ic+1,jc+1,kc]
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
                if (kc>znum(n+1)-1)
                    kc=znum(n+1)-1;
                end

                % Define normalized distances from (i,j) node to [ic,jc] node
                dx=xpos+1-jc;
                dy=ypos+1-ic;
                dz=zpos+1-kc;
        
                % Add residual from the finer level to right parts of the coarcer level 
                % for 8 nodes bounding the cell
                R(ic,jc,kc)=R(ic,jc,kc)+(1.0-dx)*(1.0-dy)*(1.0-dz)*residual(i,j,k);
                wt(ic,jc,kc)=wt(ic,jc,kc)+(1.0-dx)*(1.0-dy)*(1.0-dz);
                R(ic,jc,kc+1)=R(ic,jc,kc+1)+(1.0-dx)*(1.0-dy)*dz*residual(i,j,k);
                wt(ic,jc,kc+1)=wt(ic,jc,kc+1)+(1.0-dx)*(1.0-dy)*dz;
                R(ic+1,jc,kc)=R(ic+1,jc,kc)+(1.0-dx)*dy*(1.0-dz)*residual(i,j,k);
                wt(ic+1,jc,kc)=wt(ic+1,jc,kc)+(1.0-dx)*dy*(1.0-dz);
                R(ic+1,jc,kc+1)=R(ic+1,jc,kc+1)+(1.0-dx)*dy*dz*residual(i,j,k);
                wt(ic+1,jc,kc+1)=wt(ic+1,jc,kc+1)+(1.0-dx)*dy*dz;
                R(ic,jc+1,kc)=R(ic,jc+1,kc)+dx*(1.0-dy)*(1.0-dz)*residual(i,j,k);
                wt(ic,jc+1,kc)=wt(ic,jc+1,kc)+dx*(1.0-dy)*(1.0-dz);
                R(ic,jc+1,kc+1)=R(ic,jc+1,kc+1)+dx*(1.0-dy)*dz*residual(i,j,k);
                wt(ic,jc+1,kc+1)=wt(ic,jc+1,kc+1)+dx*(1.0-dy)*dz;
                R(ic+1,jc+1,kc)=R(ic+1,jc+1,kc)+dx*dy*(1.0-dz)*residual(i,j,k);
                wt(ic+1,jc+1,kc)=wt(ic+1,jc+1,kc)+dx*dy*(1.0-dz);
                R(ic+1,jc+1,kc+1)=R(ic+1,jc+1,kc+1)+dx*dy*dz*residual(i,j,k);
                wt(ic+1,jc+1,kc+1)=wt(ic+1,jc+1,kc+1)+dx*dy*dz;
            end
        end
    end            
end
% Recomputing right parts (R)
% for the coarcer level (k+1)
for ic=1:1:ynum(n+1);
    for jc=1:1:xnum(n+1);
        for kc=1:1:znum(n+1);
            % Check for non-zero weights and boundary nodes 
            if (wt(ic,jc,kc)==0 || bonc(ic,jc,kc)==0 || ic==1 || ic==ynum(n+1) || jc==1 || jc==xnum(n+1) || kc==1 || kc==znum(n+1))
                R(ic,jc,kc)=0;
            % Normal nodes
            else
                R(ic,jc,kc)=R(ic,jc,kc)/wt(ic,jc,kc);
            end
        end
    end            
end

