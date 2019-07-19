% Function dfidy()
% This function compose dFI/dy for specific cell boundary
% and add it to the global matrix
% Function return modified matrix
function[L]=dfidy(kp,ciA1,ciA2,ciB1,ciB2,nodx,nody,celnod,celvar,celres,pkf,L,binner,bouter,xsize,ysize)



% Upper left cell center coordinates
cxA1=(nodx(celnod(ciA1,1))+nodx(celnod(ciA1,4)))/2;
cyA1=(nody(celnod(ciA1,1))+nody(celnod(ciA1,4)))/2;
% Distance to the center of the upper left cell
cdistA1=((cxA1-xsize/2)^2+(cyA1-ysize/2)^2)^0.5;
% Lower left cell center coordinates
cxB1=(nodx(celnod(ciB1,1))+nodx(celnod(ciB1,4)))/2;
cyB1=(nody(celnod(ciB1,1))+nody(celnod(ciB1,4)))/2;
% Distance to the center of the lower left cell
cdistB1=((cxB1-xsize/2)^2+(cyB1-ysize/2)^2)^0.5;


% FIA1      
% FIB1
% Same resolution for the upper and lower cells
if(celres(ciA1)==celres(ciB1))
    if(cdistA1<binner)
        if(cdistB1<binner)
        dy=cyB1-cyA1;
        % Add  FIA1, FIB1,
        L(kp,celvar(ciA1,6))=L(kp,celvar(ciA1,6))-pkf/dy;
        L(kp,celvar(ciB1,6))=L(kp,celvar(ciB1,6))+pkf/dy;
        else
            % Lower boundary condition
            % Define vertical distance  
            % to the lower gravity potential boundary
            dy=(bouter^2-(cxA1-xsize/2)^2)^0.5+ysize/2-cyA1;
            % Add FIA1
            L(kp,celvar(ciA1,6))=L(kp,celvar(ciA1,6))-pkf/dy;
        end
    else
        % Upper boundary condition
        % Define vertical distance  
        % to the upper gravity potential boundary
        dy=(bouter^2-(cxB1-xsize/2)^2)^0.5+cyB1-ysize/2;
        % Add FIB1
        L(kp,celvar(ciB1,6))=L(kp,celvar(ciB1,6))+pkf/dy;
    end
end

% FIA1     FIA2  
%      FIB1
% Upper cells < lower cell
if(celres(ciA1)>celres(ciB1))
    % Upper right cell center coordinates
    cxA2=(nodx(celnod(ciA2,1))+nodx(celnod(ciA2,4)))/2;
    cyA2=(nody(celnod(ciA2,1))+nody(celnod(ciA2,4)))/2;
    % Distance to the center of the lower left cell
    cdistA2=((cxA2-xsize/2)^2+(cyA2-ysize/2)^2)^0.5;
    if(cdistB1<binner)
        if(cdistA1<binner && cdistA2<binner)
            dy=cyB1-cyA1;
            % Add FIB1, FIA1, FIA2
            L(kp,celvar(ciB1,6))=L(kp,celvar(ciB1,6))+pkf/dy;
            L(kp,celvar(ciA1,6))=L(kp,celvar(ciA1,6))-pkf/dy/2;
            L(kp,celvar(ciA2,6))=L(kp,celvar(ciA2,6))-pkf/dy/2;
        else
            % Upper boundary condition
            % Define vertical distance  
            % to the upper gravity potential boundary
            dy=(bouter^2-(cxB1-xsize/2)^2)^0.5+cyB1-ysize/2;
            % Add FIB1
            L(kp,celvar(ciB1,6))=L(kp,celvar(ciB1,6))+pkf/dy;
        end
    else
        % Upper left cell
        if(kp==celvar(ciA1,6))
            % Lower boundary condition
            % Define vertical distance  
            % to the lower gravity potential boundary
            dy=(bouter^2-(cxA1-xsize/2)^2)^0.5+ysize/2-cyA1;
            % Add FIA1
            L(kp,celvar(ciA1,6))=L(kp,celvar(ciA1,6))-pkf/dy;
        else
            % Upper right cell
            % Lower boundary condition
            % Define vertical distance  
            % to the lower gravity potential boundary
            dy=(bouter^2-(cxA2-xsize/2)^2)^0.5+ysize/2-cyA2;
            % Add FIA2
            L(kp,celvar(ciA2,6))=L(kp,celvar(ciA2,6))-pkf/dy;
        end
    end
end
            
%      FIA1
% FIB1     FIB2  
% Upper cell > lower cells
if(celres(ciA1)<celres(ciB1))
    % Lower right cell center coordinates
    cxB2=(nodx(celnod(ciB2,1))+nodx(celnod(ciB2,4)))/2;
    cyB2=(nody(celnod(ciB2,1))+nody(celnod(ciB2,4)))/2;
    % Distance to the center of the lower right cell
    cdistB2=((cxB2-xsize/2)^2+(cyB2-ysize/2)^2)^0.5;
    if(cdistA1<binner)
        if(cdistB1<binner && cdistB2<binner)
            dy=cyB1-cyA1;
            % Add FIA1, FIB1, FIB2
            L(kp,celvar(ciA1,6))=L(kp,celvar(ciA1,6))-pkf/dy;
            L(kp,celvar(ciB1,6))=L(kp,celvar(ciB1,6))+pkf/dy/2;
            L(kp,celvar(ciB2,6))=L(kp,celvar(ciB2,6))+pkf/dy/2;
        else
            % Lower boundary condition
            % Define vertical distance  
            % to the Right gravity potential boundary
            dy=(bouter^2-(cxA1-xsize/2)^2)^0.5+ysize/2-cyA1;
            % Add FIA1
            L(kp,celvar(ciA1,6))=L(kp,celvar(ciA1,6))-pkf/dy;
        end
    else
        % Lower left cell
        if(kp==celvar(ciB1,6))
            % Upper boundary condition
            % Define vertical distance  
            % to the upper gravity potential boundary
            dy=(bouter^2-(cxB1-xsize/2)^2)^0.5+cyB1-ysize/2;
            % Add FIB1
            L(kp,celvar(ciB1,6))=L(kp,celvar(ciB1,6))+pkf/dy;
        else
            % Lower right cell
            % Upper boundary condition
            % Define vertical distance  
            % to the upper gravity potential boundary
            dy=(bouter^2-(cxB2-xsize/2)^2)^0.5+cyB2-ysize/2;
            % Add FIB2
            L(kp,celvar(ciB2,6))=L(kp,celvar(ciB2,6))+pkf/dy;
        end
    end
end    
end

