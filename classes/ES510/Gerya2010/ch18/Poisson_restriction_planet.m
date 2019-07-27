% Function Poisson_restriction()
% This function makes restriction operation: 
% Interpolates residuals (residual) from finer (k) to coarcer (k+1) level
% using bilinear interpolation
% and produces right parts (R) for this coarcer level
%
% Resolution (xnum, ynum) and steps (xstp,ystp) for both levels
% is used for organizing the interpolation
function[R]=Poisson_restriction_planet(k,xnum,ynum,xstp,ystp,residual,bonf,bonc)

% Creating arrays for the coarser level
% Right parts
R=zeros(ynum(k+1),xnum(k+1));
% Interpolation weigts
wt=zeros(ynum(k+1),xnum(k+1));

% Interpolating residuals from finer level nodes
% Cycle of node of finer (k) level
for i=2:1:ynum(k)-1
    for j=2:1:xnum(k)-1
        % Skip nodes outside of the boundary circle
        if(bonf(i,j)>0)
            % Defining horizontal and vertical positions of current (i,j) node
            % normalized to coarcer (k+1) grid steps
            xpos=(j-1)*xstp(k)/xstp(k+1);
            ypos=(i-1)*ystp(k)/ystp(k+1);
        
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
        
            % Add residual from the finer level to right parts of the coarcer level 
            % for 4 nodes bounding the cell
            R(ic,jc)=R(ic,jc)+(1.0-dx)*(1.0-dy)*residual(i,j);
            wt(ic,jc)=wt(ic,jc)+(1.0-dx)*(1.0-dy);
            R(ic,jc+1)=R(ic,jc+1)+dx*(1.0-dy)*residual(i,j);
            wt(ic,jc+1)=wt(ic,jc+1)+dx*(1.0-dy);
            R(ic+1,jc)=R(ic+1,jc)+(1.0-dx)*dy*residual(i,j);
            wt(ic+1,jc)=wt(ic+1,jc)+(1.0-dx)*dy;
            R(ic+1,jc+1)=R(ic+1,jc+1)+dx*dy*residual(i,j);
            wt(ic+1,jc+1)=wt(ic+1,jc+1)+dx*dy;
        end
    end            
end
% Recomputing right parts (R)
% for the coarcer level (k+1)
for ic=1:1:ynum(k+1);
    for jc=1:1:xnum(k+1);
        % Check for non-zero weights and boundary nodes 
        if (wt(ic,jc)==0 || bonc(ic,jc)==0 || ic==1 || ic==ynum(k+1) || jc==1 || jc==xnum(k+1))
            R(ic,jc)=0;
        % Normal nodes
        else
            R(ic,jc)=R(ic,jc)/wt(ic,jc);
        end
    end            
end

