% Function Poisson3D_prolongation_planet()
% This function makes prolongation operation: 
% Interpolates corrections for solution (phic) from coarser (n) to finer (n-1) level
% using bilinear interpolation
% and produces updated solutions (dphi) for this finer level
%
% Resolution (xnum, ynum) and steps (xstp,ystp,zstp) for both levels
% is used for organizing the interpolation
function[dphi]=Poisson3D_prolongation_planet(n,xnum,ynum,znum,xstp,ystp,zstp,phic,bonf,bonc)

% Creating arrays for the finer level
% Solution updates
dphi=zeros(ynum(n-1),xnum(n-1),znum(n-1));
% 
% Interpolating solution updates from finer to coarser level nodes
% Cycle of node of finer (k-1) level
% No updates for boundary nodes 
for i=2:1:ynum(n-1)-1
    for j=2:1:xnum(n-1)-1
        for k=2:1:znum(n-1)-1
            % Skip nodes outside of the boundary circle
            if (bonf(i,j,k)>0)
                % Defining horizontal and vertical positions of current (i,j) node
                % normalized to coarcer (k) grid steps
                xpos=(j-1)*xstp(n-1)/xstp(n);
                ypos=(i-1)*ystp(n-1)/ystp(n);
                zpos=(k-1)*zstp(n-1)/zstp(n);
        
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
                if (jc>xnum(n)-1)
                    jc=xnum(n)-1;
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
                if (kc>znum(n)-1)
                    kc=znum(n)-1;
                end

                % Define normalized distances from (i,j) node to [ic,jc] node
                dx=xpos+1-jc;
                dy=ypos+1-ic;
                dz=zpos+1-kc;
        
                % Add solution updates from the coarser to finer level 
                % using bilinear interpolation from 8 nodes bounding the cell
                dphi(i,j,k)=dphi(i,j,k)+(1.0-dx)*(1.0-dy)*(1.0-dz)*phic(ic,jc,kc);
                dphi(i,j,k)=dphi(i,j,k)+dx*(1.0-dy)*(1.0-dz)*phic(ic,jc+1,kc);
                dphi(i,j,k)=dphi(i,j,k)+(1.0-dx)*(1.0-dz)*dy*phic(ic+1,jc,kc);
                dphi(i,j,k)=dphi(i,j,k)+dx*dy*(1.0-dz)*phic(ic+1,jc+1,kc);
                dphi(i,j,k)=dphi(i,j,k)+(1.0-dx)*(1.0-dy)*dz*phic(ic,jc,kc+1);
                dphi(i,j,k)=dphi(i,j,k)+dx*(1.0-dy)*dz*phic(ic,jc+1,kc+1);
                dphi(i,j,k)=dphi(i,j,k)+(1.0-dx)*dz*dy*phic(ic+1,jc,kc+1);
                dphi(i,j,k)=dphi(i,j,k)+dx*dy*dz*phic(ic+1,jc+1,kc+1);
            end
        end
    end
end

