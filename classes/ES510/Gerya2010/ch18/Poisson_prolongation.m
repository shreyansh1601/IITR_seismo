% Function Poisson_prolongation()
% This function makes prolongation operation: 
% Interpolates corrections for solution (phic) from coarser (k) to finer (k-1) level
% using bilinear interpolation
% and produces updated solutions (dphi) for this finer level
%
% Resolution (xnum, ynum) and steps (xstp,ystp) for both levels
% is used for organizing the interpolation
function[dphi]=Poisson_prolongation(k,xnum,ynum,xstp,ystp,phic)

% Creating arrays for the finer level
% Solution updates
dphi=zeros(ynum(k-1),xnum(k-1));
% 
% Interpolating solution updates from finer to coarser level nodes
% Cycle of node of finer (k-1) level
% No updates for boundary nodes 
for i=2:1:ynum(k-1)-1;
    for j=2:1:xnum(k-1)-1;
        % Defining horizontal and vertical positions of current (i,j) node
        % normalized to coarcer (k) grid steps
        xpos=(j-1)*xstp(k-1)/xstp(k);
        ypos=(i-1)*ystp(k-1)/ystp(k);
        
        %  [ic,jc]---------+------[ic,jc+1]
        %     |            |          |
        %     |            |          |
        %     |----------(i,j)        |
        %     |                       |
        %  [ic+1,jc]----------------[ic+1,jc+1]
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
        if (ic>ynum(k)-1)
            ic=ynum(k)-1;
        end

        % Define normalized distances from (i,j) node to [ic,jc] node
        dx=xpos+1-jc;
        dy=ypos+1-ic;
        
        % Add solution updates from the coarser to finer level 
        % using bilinear interpolation from 4 nodes bounding the cell
        dphi(i,j)=dphi(i,j)+(1.0-dx)*(1.0-dy)*phic(ic,jc);
        dphi(i,j)=dphi(i,j)+dx*(1.0-dy)*phic(ic,jc+1);
        dphi(i,j)=dphi(i,j)+(1.0-dx)*dy*phic(ic+1,jc);
        dphi(i,j)=dphi(i,j)+dx*dy*phic(ic+1,jc+1);
    end            
end

