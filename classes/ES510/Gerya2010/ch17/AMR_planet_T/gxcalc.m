% Function gxcalc()
% This function compose gx=-dFI/dx for specific cell boundary
% Function return local gx value
function[gx]=gxcalc(SP,ciL1,ciL2,ciR1,ciR2,nodx,nody,celnod,celvar,celres,binner,bouter,xsize,ysize)



% Upper left cell center coordinates
cxL1=(nodx(celnod(ciL1,1))+nodx(celnod(ciL1,4)))/2;
cyL1=(nody(celnod(ciL1,1))+nody(celnod(ciL1,4)))/2;
% Distance to the center of the upper left cell
cdistL1=((cxL1-xsize/2)^2+(cyL1-ysize/2)^2)^0.5;
% Upper right cell center coordinates
cxR1=(nodx(celnod(ciR1,1))+nodx(celnod(ciR1,4)))/2;
cyR1=(nody(celnod(ciR1,1))+nody(celnod(ciR1,4)))/2;
% Distance to the center of the upper right cell
cdistR1=((cxR1-xsize/2)^2+(cyR1-ysize/2)^2)^0.5;


% FIL1      FIR1
% Same resolution for the left and right cells
if(celres(ciL1)==celres(ciR1))
    if(cdistL1<binner)
        if(cdistR1<binner)
        dx=cxR1-cxL1;
        % Add  FIL1, FIR1,
        gx=(SP(celvar(ciL1,6))-SG(celvar(ciR1,6)))/dx;
        else
            % Right boundary condition
            % Define horizontal distance  
            % to the right gravity potential boundary
            dx=(bouter^2-(cyL1-ysize/2)^2)^0.5+xsize/2-cxL1;
            % Add FIL1
            gx=SP(celvar(ciL1,6))/dx;
        end
    else
        % Left boundary condition
        % Define horizontal distance  
        % to the Left gravity potential boundary
        dx=(bouter^2-(cyR1-ysize/2)^2)^0.5+cxR1-xsize/2;
        % Add FIR1
        gx=-SG(celvar(ciR1,6))/dx;
    end
end

% FIL1      
%           FIR1
% FIL2
% Left cells < right cell
if(celres(ciL1)>celres(ciR1))
    % Lower left cell center coordinates
    cxL2=(nodx(celnod(ciL2,1))+nodx(celnod(ciL2,4)))/2;
    cyL2=(nody(celnod(ciL2,1))+nody(celnod(ciL2,4)))/2;
    % Distance to the center of the lower left cell
    cdistL2=((cxL2-xsize/2)^2+(cyL2-ysize/2)^2)^0.5;
    if(cdistR1<binner)
        if(cdistL1<binner && cdistL2<binner)
            dx=cxR1-cxL1;
            % Add FIR1, FIL1, FIL2
            gx=((SG(celvar(ciL1,6))+SG(celvar(ciL1,6)))/-SP(celvar(ciR1,6)))/dx;
        else
            % Left boundary condition
            % Define horizontal distance  
            % to the Left gravity potential boundary
            dx=(bouter^2-(cyR1-ysize/2)^2)^0.5+cxR1-xsize/2;
            % Add FIR1
            gx=-SG(celvar(ciR1,6))/dx;
        end
    else
        % Upper left cell
        if(kp==celvar(ciL1,6))
            % Right boundary condition
            % Define horizontal distance  
            % to the right gravity potential boundary
            dx=(bouter^2-(cyL1-ysize/2)^2)^0.5+xsize/2-cxL1;
            % Add FIL1
            L(kp,celvar(ciL1,6))=L(kp,celvar(ciL1,6))-pkf/dx;
        else
            % Lower left cell
            % Right boundary condition
            % Define horizontal distance  
            % to the right gravity potential boundary
            dx=(bouter^2-(cyL2-ysize/2)^2)^0.5+xsize/2-cxL2;
            % Add FIL2
            L(kp,celvar(ciL2,6))=L(kp,celvar(ciL2,6))-pkf/dx;
        end
    end
end
            
%           FIR1      
% FIL1
%           FIR2
% Left cell > right cells
if(celres(ciL1)<celres(ciR1))
    % Lower right cell center coordinates
    cxR2=(nodx(celnod(ciR2,1))+nodx(celnod(ciR2,4)))/2;
    cyR2=(nody(celnod(ciR2,1))+nody(celnod(ciR2,4)))/2;
    % Distance to the center of the lower right cell
    cdistR2=((cxR2-xsize/2)^2+(cyR2-ysize/2)^2)^0.5;
    if(cdistL1<binner)
        if(cdistR1<binner && cdistR2<binner)
            dx=cxR1-cxL1;
            % Add FIL1, FIR1, FIR2
            L(kp,celvar(ciL1,6))=L(kp,celvar(ciL1,6))-pkf/dx;
            L(kp,celvar(ciR1,6))=L(kp,celvar(ciR1,6))+pkf/dx/2;
            L(kp,celvar(ciR2,6))=L(kp,celvar(ciR2,6))+pkf/dx/2;
        else
            % Right boundary condition
            % Define horizontal distance  
            % to the Right gravity potential boundary
            dx=(bouter^2-(cyL1-ysize/2)^2)^0.5+xsize/2-cxL1;
            % Add FIL1
            L(kp,celvar(ciL1,6))=L(kp,celvar(ciL1,6))-pkf/dx;
        end
    else
        % Upper right cell
        if(kp==celvar(ciR1,6))
            % Left boundary condition
            % Define horizontal distance  
            % to the left gravity potential boundary
            dx=(bouter^2-(cyR1-ysize/2)^2)^0.5+cxR1-xsize/2;
            % Add FIR1
            L(kp,celvar(ciR1,6))=L(kp,celvar(ciR1,6))+pkf/dx;
        else
            % Lower right cell
            % Left boundary condition
            % Define horizontal distance  
            % to the left gravity potential boundary
            dx=(bouter^2-(cyR2-ysize/2)^2)^0.5+cxR2-xsize/2;
            % Add FIR2
            L(kp,celvar(ciR2,6))=L(kp,celvar(ciR2,6))+pkf/dx;
        end
    end
end    
end

