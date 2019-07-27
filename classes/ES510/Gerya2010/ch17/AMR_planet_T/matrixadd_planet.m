% Function matrixad_planet()
% This function compose matrix and right parts
% for Poisson and then for continuity and Stokes equations
% Function return solutions for Poisson SP() and Stokes S()
function[S,SP,SG]=matrixadd_planet(timestp,varnum,varnump,celnum,celdot,celnod,celvar,celres,celeta,celrho,nodx,nody,nodcel,nodeta,nodrho,bupper,blower,bleft,bright,binner,bouter,xsize,ysize,pscale)
% timestp=0;

% Define global matrices for Poisson
L=sparse(varnump,varnump);
R=zeros(varnump,1);

% Composing equations
for ci=1:1:celnum
    if(celdot(ci,1)==0)
        %
        %       |ciA1 |ciA2 | 
        %   ----+-----+-----+----
        %  ciL1 |  ci       | ciR1
        %   ----+     FI    +---
        %  ciL2 |           | ciR2
        %   ----+-----+-----+----
        %       |ciB1 |ciB2 |
        % Indexes of cells
        % Left
        ciL1=nodcel(celnod(ci,1),2); 
        ciL2=nodcel(celnod(ci,2),1);
        % Right
        ciR1=nodcel(celnod(ci,3),4);
        ciR2=nodcel(celnod(ci,4),3);
        % Above
        ciA1=nodcel(celnod(ci,1),3); 
        ciA2=nodcel(celnod(ci,3),1); 
        % Below
        ciB1=nodcel(celnod(ci,2),4);
        ciB2=nodcel(celnod(ci,4),2);
        
        % Equation for FI------------------------
        % d2G/dx2+d2G/dy2=4*KF*PI*RHO
        % Index for FI
        kp=celvar(ci,6);
        % Cell sizes
        dx=(nodx(celnod(ci,3))-nodx(celnod(ci,1))); 
        dy=(nody(celnod(ci,2))-nody(celnod(ci,1))); 
         % Cell center coordinates
        cx=(nodx(celnod(ci,1))+nodx(celnod(ci,4)))/2;
        cy=(nody(celnod(ci,1))+nody(celnod(ci,4)))/2;
        % Distance to the center
        cdist=((cx-xsize/2)^2+(cy-ysize/2)^2)^0.5;
        
        % Boundary condition
        if(cdist>=binner)
            L(kp,kp)=1;
            R(kp)=0;
        else
            % Poisson equation
            %                 FIA
            %               dFI/dy1
            %  FIL  dFI/dx1    FI   dFI/dx2  FIR
            %               dFI/dy2
            %                 FIB
            % Compose matrix
            % d2G/dx2+d2G/dy2=4*KF*GAMMA*PI*RHO
            % -dgx/dx-dgy/dy=4*KF*GAMMA*PI*RHO
            %
            % Add matrix for dFI/dx1
            if(ciL1>0)
                [L]=dfidx(kp,ciL1,ciL2,ci,ciB1,nodx,nody,celnod,celvar,celres,-1/dx,L,binner,bouter,xsize,ysize);
            else
                [L]=dfidx(kp,ciL2,ciL2,ciA1,ci,nodx,nody,celnod,celvar,celres,-1/dx,L,binner,bouter,xsize,ysize);
            end
            %
            % Add matrix for dFI/dx2
            if(ciR1>0)
                [L]=dfidx(kp,ci,ciB2,ciR1,ciR2,nodx,nody,celnod,celvar,celres,1/dx,L,binner,bouter,xsize,ysize);
            else
                [L]=dfidx(kp,ciA2,ci,ciR2,ciR2,nodx,nody,celnod,celvar,celres,1/dx,L,binner,bouter,xsize,ysize);
            end
            %
            % Add matrix for dFI/dy1
            if(ciA1>0)
                [L]=dfidy(kp,ciA1,ciA2,ci,ciR1,nodx,nody,celnod,celvar,celres,-1/dy,L,binner,bouter,xsize,ysize);
            else
                [L]=dfidy(kp,ciA2,ciA2,ciL1,ci,nodx,nody,celnod,celvar,celres,-1/dy,L,binner,bouter,xsize,ysize);
            end
            %
            % Add matrix for dFI/dy2
            if(ciB1>0)
                [L]=dfidy(kp,ci,ciR2,ciB1,ciB2,nodx,nody,celnod,celvar,celres,1/dy,L,binner,bouter,xsize,ysize);
            else
                [L]=dfidy(kp,ciL2,ci,ciB2,ciB2,nodx,nody,celnod,celvar,celres,1/dy,L,binner,bouter,xsize,ysize);
            end
            %
            % Right part
            % GK=6.672e-11 Gravity constant, N*m^2/kg^2 
            R(kp)=4*2/3*6.672e-11*pi*celrho(ci);
        end
    end
end


% % Showing matrix structure
% figure(4)
% spy(L);
% 
% Solving the matrix
SP=L\R;

% Define global matrices for Stokes
L=sparse(varnum,varnum);
R=zeros(varnum,1);
% Processed Y/N index
varyn=zeros(varnum,1);
% Array for gravity acceleration components
SG=zeros(varnum,1);

% timestp=0;

% Composing equations
for ci=1:1:celnum
    if(celdot(ci,1)==0)
        %
        %       |ciA1 |ciA2 | 
        %   ----+----vy2----+----
        %  ciL1 |  ci       | ciR1
        %   ---vx1   P5    vx4---
        %  ciL2 |           | ciR2
        %   ----+----vy3----+----
        %       |ciB1 |ciB2 |
        % Indexes of cells
        % Left
        ciL1=nodcel(celnod(ci,1),2); 
        ciL2=nodcel(celnod(ci,2),1);
        % Right
        ciR1=nodcel(celnod(ci,3),4);
        ciR2=nodcel(celnod(ci,4),3);
        % Above
        ciA1=nodcel(celnod(ci,1),3); 
        ciA2=nodcel(celnod(ci,3),1); 
        % Below
        ciB1=nodcel(celnod(ci,2),4);
        ciB2=nodcel(celnod(ci,4),2);
        
        % Equation for P5------------------------
        % Index for P5
        kp=celvar(ci,5);
        % Boundary condition
        if(nodx(celnod(ci,1))==0 && nody(celnod(ci,1))==0)
            L(kp,kp)=1*pscale;
            R(kp)=0;
        else
            % Continuity equation
            %     vy2
            % vx1  P5   vx4
            %     vy3
            kfx=1/(nodx(celnod(ci,3))-nodx(celnod(ci,1)));
            kfy=1/(nody(celnod(ci,2))-nody(celnod(ci,1)));
            L(kp,celvar(ci,1))=-kfx; % vx1
            L(kp,celvar(ci,2))=-kfy; % vy2
            L(kp,celvar(ci,3))=+kfy; % vy3
            L(kp,celvar(ci,4))=+kfx; % vx4
            R(kp)=0;
        end
        % Mark processed equation
        varyn(kp)=1;
        % End Equation for P5------------------------
        
        
        % Equation for vx1------------------------
        % Index for vx1
        kx=celvar(ci,1);
        % Check processed variable y/n
        if(varyn(kx)==0)
            % Boundary condition
            if(nodx(celnod(ci,1))==0)
                L(kx,kx)=1;
                R(kx)=0;
                % Mark processed equation
                varyn(kx)=1;
            else
                % Solving x-Stokes
                % dSIGxx/dx-dP/dx + dSIGxy/dy - dt/2*gx*(vx*dRHO/dx+vy*dRHO/dy)= -RHO*gx
                % Cell/Subcells -> Cell   
                if(ciL1>0 && celres(ciL1)>=celres(ci))
                    % Compose x-Stokes equation 
                    % dSIGxx/dx-dP/dx + dSIGxy/dy = -RHO*gx
                    %            SIGxy1
                    % SIG'xx1,P1    vx1    SIG'xx2,P2    vx4
                    %            SIGxy2
                    % Compose matrix
                    % dSIGxx/dx-dP/dx
                    % x-distance between SIGxx-nodes
                    dx=(nodx(celnod(ci,3))-nodx(celnod(ciL1,1)))/2; 
                    % Add SIGxx2=SIG'xx2-P2
                    [L]=sigmaxx(kx,ci,nodx,celnod,celvar,celeta,1/dx,pscale/dx,L);
                    % Add SIGxx1=SIG'xx1-P1
                    if(celres(ciL1)==celres(ci))
                        % No change in resolution
                        [L]=sigmaxx(kx,ciL1,nodx,celnod,celvar,celeta,-1/dx,-pscale/dx,L);
                        % Compute gx=-dFI/dx
%                         if(SP(celvar(ci,6))==0 || SP(celvar(ciL1,6))==0) 
%                             SG(celvar(ci,1))=0;
%                         else
                            SG(celvar(ci,1))=-(SP(celvar(ci,6))-SP(celvar(ciL1,6)))/dx;
%                         end 
                        % Add -gx*vx*dt/2*dRHO/dx
                        L(kx,kx)=L(kx,kx)-SG(celvar(ci,1))*timestp/2*(celrho(ci)-celrho(ciL1))/dx;
                    else
                        % Increase in resolution
                        % Upper left cell
                        [L]=sigmaxx(kx,ciL1,nodx,celnod,celvar,celeta,-0.5/dx,-0.5*pscale/dx,L);                        % P1
                        % Lower left cell
                        [L]=sigmaxx(kx,ciL2,nodx,celnod,celvar,celeta,-0.5/dx,-0.5*pscale/dx,L);                        % P1
                        % Compute gx=-dFI/dx
%                         if(SP(celvar(ci,6))==0 || SP(celvar(ciL1,6))==0 || SP(celvar(ciL2,6))==0) 
%                             SG(celvar(ci,1))=0;
%                         else
                            SG(celvar(ci,1))=-(SP(celvar(ci,6))-(SP(celvar(ciL1,6))+SP(celvar(ciL2,6)))/2)/dx;
%                         end
                        % Add -gx*vx*dt/2*dRHO/dx
                        L(kx,kx)=L(kx,kx)-SG(celvar(ci,1))*timestp/2*(celrho(ci)-(celrho(ciL1)+celrho(ciL2))/2)/dx;
                    end
                    % dSIGxy/dy
                    % y-distance between SIGxy-nodes
                    dy=nody(celnod(ci,2))-nody(celnod(ci,1)); 
                    % Add -gx*vy*dt/2*dRHO/dy
                    ival=-SG(celvar(ci,1))*timestp/8*(nodrho(celnod(ci,2))-nodrho(celnod(ci,1)))/dy;
                    L(kx,celvar(ci,2))=L(kx,celvar(ci,2))+ival;
                    L(kx,celvar(ci,3))=L(kx,celvar(ci,3))+ival;
                    L(kx,celvar(ciL1,2))=L(kx,celvar(ciL1,2))+ival;
                    L(kx,celvar(ciL2,3))=L(kx,celvar(ciL2,3))+ival;
                    % Add SIGxy1
                    if(ciA1==0 && ciA2>0)
                        % vx between two subcells
                        %   |             |
                        % -ni1-----+-----ni2-
                        %   |      |      |       
                        %         vx1
                        % Left top node
                        [L]=sigmaxy(kx,celnod(ciL1,1),nodx,nody,nodcel,celnod,celvar,celres,nodeta,-0.5/dy,L,bupper,blower,bleft,bright,xsize,ysize);
                        % Right top node
                        [L]=sigmaxy(kx,celnod(ci,3),nodx,nody,nodcel,celnod,celvar,celres,nodeta,-0.5/dy,L,bupper,blower,bleft,bright,xsize,ysize);
                    else
                        % vx in standard position
                        %          |
                        %   ------ni1------
                        %          |             
                        %         vx1
                        % Top node
                        [L]=sigmaxy(kx,celnod(ci,1),nodx,nody,nodcel,celnod,celvar,celres,nodeta,-1/dy,L,bupper,blower,bleft,bright,xsize,ysize);
                    end
                    % Add SIGxy2
                    if(ciB1==0 && ciB2>0)
                        % vx between two subcells
                        %         vx1
                        %   |      |      |       
                        % -ni1-----+-----ni2-
                        %   |             |
                        % Left bottom node
                        [L]=sigmaxy(kx,celnod(ciL2,2),nodx,nody,nodcel,celnod,celvar,celres,nodeta,0.5/dy,L,bupper,blower,bleft,bright,xsize,ysize);
                        % Right bottom node
                        [L]=sigmaxy(kx,celnod(ci,4),nodx,nody,nodcel,celnod,celvar,celres,nodeta,0.5/dy,L,bupper,blower,bleft,bright,xsize,ysize);
                    else
                        % vx in standard position
                        %         vx1
                        %          |
                        %   ------ni1------
                        %          |             
                        % Top node
                        [L]=sigmaxy(kx,celnod(ci,2),nodx,nody,nodcel,celnod,celvar,celres,nodeta,1/dy,L,bupper,blower,bleft,bright,xsize,ysize);
                    end
                    % Right Part -RHO*gx
                    R(kx)=-SG(celvar(ci,1))*(nodrho(celnod(ci,1))+nodrho(celnod(ci,2)))/2;
                    % Mark processed equation
                    varyn(kx)=1;
                end
            end
        end
        % End Equation for vx1------------------------
        
        
        % Equation for vx4------------------------
        % Index for vx4
        kx=celvar(ci,4);
        % Check processed variable y/n
        if(varyn(kx)==0)
            % Boundary condition
            if(nodx(celnod(ci,3))==xsize)
                L(kx,kx)=1;
                R(kx)=0;
                % Mark processed equation
                varyn(kx)=1;
           else
                % Solving x-Stokes
                % dSIGxx/dx-dP/dx + dSIGxy/dy - dt/2*gx*(vx*dRHO/dx+vy*dRHO/dy)= -RHO*gx
                %  Cell->Cell/Subcells    
                if(ciR1>0 && celres(ciR1)>=celres(ci))
                    % Compose x-Stokes equation 
                    % dSIGxx/dx-dP/dx + dSIGxy/dy = -RHO*gx
                    %                    SIGxy1
                    % vx1   SIG'xx1,P1    vx4    SIG'xx2,P2    
                    %                    SIGxy2
                    % Compose matrix
                    % dSIGxx/dx-dP/dx
                    % x-distance between SIGxx-nodes
                    dx=(nodx(celnod(ciR1,3))-nodx(celnod(ci,1)))/2; 
                    % Add SIGxx1=SIG'xx1-P1
                    [L]=sigmaxx(kx,ci,nodx,celnod,celvar,celeta,-1/dx,-pscale/dx,L);
                    % Add SIGxx2=SIG'xx2-P2
                    if(celres(ciR1)==celres(ci))
                        % No change in resolution
                        [L]=sigmaxx(kx,ciR1,nodx,celnod,celvar,celeta,1/dx,pscale/dx,L);
                        % Compute gx=-dFI/dx
%                         if(SP(celvar(ci,6))==0 || SP(celvar(ciR1,6))==0) 
%                             SG(celvar(ci,4))=0;
%                         else
                            SG(celvar(ci,4))=-(SP(celvar(ciR1,6))-SP(celvar(ci,6)))/dx;
%                         end
                        % Add -gx*vx*dt/2*dRHO/dx
                        L(kx,kx)=L(kx,kx)-SG(celvar(ci,4))*timestp/2*(celrho(ciR1)-celrho(ci))/dx;
                    else
                        % Increase in resolution
                        % Upper right cell
                        [L]=sigmaxx(kx,ciR1,nodx,celnod,celvar,celeta,0.5/dx,0.5*pscale/dx,L);                        % P1
                        % Lower right cell
                        [L]=sigmaxx(kx,ciR2,nodx,celnod,celvar,celeta,0.5/dx,0.5*pscale/dx,L);                        % P1
                        % Compute gx=-dFI/dx
%                         if(SP(celvar(ci,6))==0 || SP(celvar(ciR1,6))==0 || SP(celvar(ciR2,6))==0) 
%                             SG(celvar(ci,4))=0;
%                         else
                            SG(celvar(ci,4))=-((SP(celvar(ciR1,6))+SP(celvar(ciR2,6)))/2-SP(celvar(ci,6)))/dx;
%                         end
                        % Add -gx*vx*dt/2*dRHO/dx
                        L(kx,kx)=L(kx,kx)-SG(celvar(ci,4))*timestp/2*((celrho(ciR1)+celrho(ciR2))/2-celrho(ci))/dx;
                    end
                    % dSIGxy/dy
                    % y-distance between SIGxy-nodes
                    dy=nody(celnod(ci,4))-nody(celnod(ci,3)); 
                    % Add -gx*vy*dt/2*dRHO/dy
                    ival=-SG(celvar(ci,4))*timestp/8*(nodrho(celnod(ci,4))-nodrho(celnod(ci,3)))/dy;
                    L(kx,celvar(ci,2))=L(kx,celvar(ci,2))+ival;
                    L(kx,celvar(ci,3))=L(kx,celvar(ci,3))+ival;
                    L(kx,celvar(ciR1,2))=L(kx,celvar(ciR1,2))+ival;
                    L(kx,celvar(ciR2,3))=L(kx,celvar(ciR2,3))+ival;
                    % Add SIGxy1
                    if(ciA1>0 && ciA2==0)
                        % vx between two subcells
                        %   |     ciA1    |
                        % -ni1-----+-----ni2-
                        %   |  ci  |      |       
                        %         vx4
                        % Left top node
                        [L]=sigmaxy(kx,celnod(ci,1),nodx,nody,nodcel,celnod,celvar,celres,nodeta,-0.5/dy,L,bupper,blower,bleft,bright,xsize,ysize);
                        % Right top node
                        [L]=sigmaxy(kx,celnod(ciR1,3),nodx,nody,nodcel,celnod,celvar,celres,nodeta,-0.5/dy,L,bupper,blower,bleft,bright,xsize,ysize);
                    else
                        % vx in standard position
                        %          |
                        %   ------ni1------
                        %          |             
                        %         vx4
                        % Top node
                        [L]=sigmaxy(kx,celnod(ci,3),nodx,nody,nodcel,celnod,celvar,celres,nodeta,-1/dy,L,bupper,blower,bleft,bright,xsize,ysize);
                    end
                    % Add SIGxy2
                    if(ciB1>0 && ciB2==0)
                        % vx between two subcells
                        %         vx4
                        %   |  ci  |      |       
                        % -ni1-----+-----ni2-
                        %   |     ciB1    |
                        % Left bottom node
                        [L]=sigmaxy(kx,celnod(ci,2),nodx,nody,nodcel,celnod,celvar,celres,nodeta,0.5/dy,L,bupper,blower,bleft,bright,xsize,ysize);
                        % Right bottom node
                        [L]=sigmaxy(kx,celnod(ciR2,4),nodx,nody,nodcel,celnod,celvar,celres,nodeta,0.5/dy,L,bupper,blower,bleft,bright,xsize,ysize);
                    else
                        % vx in standard position
                        %         vx1
                        %          |
                        %   ------ni1------
                        %          |             
                        % Bottom node
                        [L]=sigmaxy(kx,celnod(ci,4),nodx,nody,nodcel,celnod,celvar,celres,nodeta,1/dy,L,bupper,blower,bleft,bright,xsize,ysize);
                    end
                    % Right Part
                    R(kx)=-SG(celvar(ci,4))*(nodrho(celnod(ci,3))+nodrho(celnod(ci,4)))/2;
                    % Mark processed equation
                    varyn(kx)=1;
                end
            end
        end
        % End Equation for vx4------------------------
        
        
        % Equation for vy2------------------------
        % Index for vy2
        ky=celvar(ci,2);
        % Check processed variable y/n
        if(varyn(ky)==0)
            % Boundary condition
            if(nody(celnod(ci,1))==0)
                L(ky,ky)=1;
                R(ky)=0;
                % Mark processed equation
                varyn(ky)=1;
            else
                % Solving y-Stokes
                % dSIGyy/dy-dP/dy + dSIGxy/dx - dt/2*gy*(vx*dRHO/dx+vy*dRHO/dy) = -RHO*gy
                % Cell/Subcells -> Cell   
                if(ciA1>0 && celres(ciA1)>=celres(ci))
                    % Compose y-Stokes equation 
                    % dSIGyy/dy-dP/dy + dSIGxy/dx = -RHO*gy
                    %           SIGyy1,P1
                    %   SIGxy1     vy2    SIG'xy2
                    %           SIGyy2,P2
                    %              vy3
                    % Compose matrix
                    % dSIGyy/dy-dP/dy
                    % y-distance between SIGyy-nodes
                    dy=(nody(celnod(ci,2))-nody(celnod(ciA1,1)))/2; 
                    % Add SIGyy2=SIG'yy2-P2
                    [L]=sigmayy(ky,ci,nody,celnod,celvar,celeta,1/dy,pscale/dy,L);
                    % Add SIGyy1=SIG'yy1-P1
                    if(celres(ciA1)==celres(ci))
                        % No change in resolution
                        [L]=sigmayy(ky,ciA1,nody,celnod,celvar,celeta,-1/dy,-pscale/dy,L);
                        % Compute gy=-dFI/dy
%                         if(SP(celvar(ci,6))==0 || SP(celvar(ciA1,6))==0) 
%                             SG(celvar(ci,2))=0;
%                         else
                            SG(celvar(ci,2))=-(SP(celvar(ci,6))-SP(celvar(ciA1,6)))/dy;
%                         end
                        % Add -gy*vy*dt/2*dRHO/dy
                        L(ky,ky)=L(ky,ky)-SG(celvar(ci,2))*timestp/2*(celrho(ci)-celrho(ciA1))/dy;
                    else
                        % Increase in resolution
                        % Upper left cell
                        [L]=sigmayy(ky,ciA1,nody,celnod,celvar,celeta,-0.5/dy,-0.5*pscale/dy,L);                        % P1
                        % Upper right cell
                        [L]=sigmayy(ky,ciA2,nody,celnod,celvar,celeta,-0.5/dy,-0.5*pscale/dy,L);                        % P1
                        % Compute gy=-dFI/dy
%                         if(SP(celvar(ci,6))==0 || SP(celvar(ciA1,6))==0 || SP(celvar(ciA2,6))==0) 
%                             SG(celvar(ci,2))=0;
%                         else
                            SG(celvar(ci,2))=-(SP(celvar(ci,6))-(SP(celvar(ciA1,6))+SP(celvar(ciA2,6)))/2)/dy;
%                         end
                        % Add -gy*vy*dt/2*dRHO/dy
                        L(ky,ky)=L(ky,ky)-SG(celvar(ci,2))*timestp/2*(celrho(ci)-(celrho(ciA1)+celrho(ciA2))/2)/dy;
                   end
                    % dSIGxy/dx
                    % x-distance between SIGxy-nodes
                    dx=nodx(celnod(ci,3))-nodx(celnod(ci,1)); 
                     % Add -gy*vx*dt/2*dRHO/dx
                    ival=-SG(celvar(ci,2))*timestp/8*(nodrho(celnod(ci,3))-nodrho(celnod(ci,1)))/dx;
                    L(ky,celvar(ci,1))=L(ky,celvar(ci,1))+ival;
                    L(ky,celvar(ci,4))=L(ky,celvar(ci,4))+ival;
                    L(ky,celvar(ciA1,1))=L(ky,celvar(ciA1,1))+ival;
                    L(ky,celvar(ciA2,4))=L(ky,celvar(ciA2,4))+ival;
                   % Add SIGxy1
                    if(ciL1==0 && ciL2>0)
                        % vy between two subcells
                        %       |
                        %------ni1------
                        %       | ciA1
                        %  ciL2 +--vy2
                        %       | ci
                        %------ni2------
                        %       |
                        % Upper left node
                        [L]=sigmaxy(ky,celnod(ciA1,1),nodx,nody,nodcel,celnod,celvar,celres,nodeta,-0.5/dx,L,bupper,blower,bleft,bright,xsize,ysize);
                        % Lower left node
                        [L]=sigmaxy(ky,celnod(ci,2),nodx,nody,nodcel,celnod,celvar,celres,nodeta,-0.5/dx,L,bupper,blower,bleft,bright,xsize,ysize);
                    else
                        % vy in standard position
                        %          |
                        %   ------ni1---vy2
                        %          |  ci           
                        % Left node
                        [L]=sigmaxy(ky,celnod(ci,1),nodx,nody,nodcel,celnod,celvar,celres,nodeta,-1/dx,L,bupper,blower,bleft,bright,xsize,ysize);
                    end
                    % Add SIGxy2
                    if(ciR1==0 && ciR2>0)
                        % vy between two subcells
                        %         |
                        %  ------ni1------
                        %    ciA2 | 
                        %  --vy2--+ ciR2
                        %     ci  | 
                        %  ------ni2------
                        %         |
                        % Left bottom node
                        [L]=sigmaxy(ky,celnod(ciA2,3),nodx,nody,nodcel,celnod,celvar,celres,nodeta,0.5/dx,L,bupper,blower,bleft,bright,xsize,ysize);
                        % Right bottom node
                        [L]=sigmaxy(ky,celnod(ci,4),nodx,nody,nodcel,celnod,celvar,celres,nodeta,0.5/dx,L,bupper,blower,bleft,bright,xsize,ysize);
                    else
                        % vy in standard position
                        %          |
                        %   vy2---ni1---
                        %     ci   |            
                        % Top node
                        [L]=sigmaxy(ky,celnod(ci,3),nodx,nody,nodcel,celnod,celvar,celres,nodeta,1/dx,L,bupper,blower,bleft,bright,xsize,ysize);
                    end
                    % Right Part
                    R(ky)=-SG(celvar(ci,2))*(nodrho(celnod(ci,1))+nodrho(celnod(ci,3)))/2;
                    % Mark processed equation
                    varyn(ky)=1;
                end
            end
        end
        % End Equation for vy2------------------------
                
        
        
        % Equation for vy3------------------------
        % Index for vy3
        ky=celvar(ci,3);
        % Check processed variable y/n
        if(varyn(ky)==0)
            % Boundary condition
            if(nody(celnod(ci,4))==ysize)
                L(ky,ky)=1;
                R(ky)=0;
                % Mark processed equation
                varyn(ky)=1;
            else
                % Solving y-Stokes
                % dSIGyy/dy-dP/dy + dSIGxy/dx - dt/2*gy*(vx*dRHO/dx+vy*dRHO/dy) = -RHO*gy
                % Cell/Subcells -> Cell   
                if(ciB1>0 && celres(ciB1)>=celres(ci))
                    % Compose y-Stokes equation 
                    % dSIGyy/dy-dP/dy + dSIGxy/dx = -RHO*gy
                    %              vy2
                    %           SIGyy1,P1
                    %   SIGxy1     vy3    SIG'xy2
                    %           SIGyy2,P2
                    % Compose matrix
                    % dSIGyy/dy-dP/dy
                    % y-distance between SIGyy-nodes
                    dy=(nody(celnod(ciB1,2))-nody(celnod(ci,1)))/2; 
                    % Add SIGyy1=SIG'yy1-P1
                    [L]=sigmayy(ky,ci,nody,celnod,celvar,celeta,-1/dy,-pscale/dy,L);
                    % Add SIGyy2=SIG'yy2-P2
                    if(celres(ciB1)==celres(ci))
                        % No change in resolution
                        [L]=sigmayy(ky,ciB1,nody,celnod,celvar,celeta,1/dy,pscale/dy,L);
                        % Compute gy=-dFI/dy
%                         if(SP(celvar(ci,6))==0 || SP(celvar(ciB1,6))==0) 
%                             SG(celvar(ci,3))=0;
%                         else
                            SG(celvar(ci,3))=-(SP(celvar(ciB1,6))-SP(celvar(ci,6)))/dy;
%                         end
                        % Add -gy*vy*dt/2*dRHO/dy
                        L(ky,ky)=L(ky,ky)-SG(celvar(ci,3))*timestp/2*(celrho(ciB1)-celrho(ci))/dy;
                   else
                        % Increase in resolution
                        % Upper left cell
                        [L]=sigmayy(ky,ciB1,nody,celnod,celvar,celeta,0.5/dy,0.5*pscale/dy,L);                        % P1
                        % Upper right cell
                        [L]=sigmayy(ky,ciB2,nody,celnod,celvar,celeta,0.5/dy,0.5*pscale/dy,L);                        % P1
                        % Compute gy=-dFI/dy
%                         if(SP(celvar(ci,6))==0 || SP(celvar(ciB1,6))==0 || SP(celvar(ciB2,6))==0) 
%                             SG(celvar(ci,3))=0;
%                         else
                            SG(celvar(ci,3))=-((SP(celvar(ciB1,6))+SP(celvar(ciB2,6)))/2-SP(celvar(ci,6)))/dy;
%                         end
                        % Add -gy*vy*dt/2*dRHO/dy
                        L(ky,ky)=L(ky,ky)-SG(celvar(ci,3))*timestp/2*((celrho(ciB1)+celrho(ciB2))/2-celrho(ci))/dy;
                    end
                    % dSIGxy/dx
                    % x-distance between SIGxy-nodes
                    dx=nodx(celnod(ci,4))-nodx(celnod(ci,2)); 
                     % Add -gy*vx*dt/2*dRHO/dx
                    ival=-SG(celvar(ci,3))*timestp/8*(nodrho(celnod(ci,4))-nodrho(celnod(ci,2)))/dx;
                    L(ky,celvar(ci,1))=L(ky,celvar(ci,1))+ival;
                    L(ky,celvar(ci,4))=L(ky,celvar(ci,4))+ival;
                    L(ky,celvar(ciB1,1))=L(ky,celvar(ciB1,1))+ival;
                    L(ky,celvar(ciB2,4))=L(ky,celvar(ciB2,4))+ival;
                    % Add SIGxy1
                    if(ciL2==0 && ciL1>0)
                        % vy between two subcells
                        %       |
                        %------ni1------
                        %       | ci
                        %  ciL1 +--vy3
                        %       | ciB1
                        %------ni2------
                        %       |
                        % Upper left node
                        [L]=sigmaxy(ky,celnod(ci,1),nodx,nody,nodcel,celnod,celvar,celres,nodeta,-0.5/dx,L,bupper,blower,bleft,bright,xsize,ysize);
                        % Lower left node
                        [L]=sigmaxy(ky,celnod(ciB1,2),nodx,nody,nodcel,celnod,celvar,celres,nodeta,-0.5/dx,L,bupper,blower,bleft,bright,xsize,ysize);
                    else
                        % vy in standard position
                        %          |   ci 
                        %   ------ni1---vy3
                        %          |            
                        % Left node
                        [L]=sigmaxy(ky,celnod(ci,2),nodx,nody,nodcel,celnod,celvar,celres,nodeta,-1/dx,L,bupper,blower,bleft,bright,xsize,ysize);
                    end
                    % Add SIGxy2
                    if(ciR2==0 && ciR1>0)
                        % vy between two subcells
                        %         |
                        %  ------ni1------
                        %    ci   | 
                        %  --vy3--+ ciR1
                        %    ciB2 | 
                        %  ------ni2------
                        %         |
                        % Right top node
                        [L]=sigmaxy(ky,celnod(ci,3),nodx,nody,nodcel,celnod,celvar,celres,nodeta,0.5/dx,L,bupper,blower,bleft,bright,xsize,ysize);
                        % Right bottom node
                        [L]=sigmaxy(ky,celnod(ciB2,4),nodx,nody,nodcel,celnod,celvar,celres,nodeta,0.5/dx,L,bupper,blower,bleft,bright,xsize,ysize);
                    else
                        % vy in standard position
                        %          |
                        %   vy3---ni1---
                        %     ci   |            
                        % Top node
                        [L]=sigmaxy(ky,celnod(ci,4),nodx,nody,nodcel,celnod,celvar,celres,nodeta,1/dx,L,bupper,blower,bleft,bright,xsize,ysize);
                    end
                    % Right Part
                    R(ky)=-SG(celvar(ci,3))*(nodrho(celnod(ci,2))+nodrho(celnod(ci,4)))/2;
                    % Mark processed equation
                    varyn(ky)=1;
                end
            end
        end
        % End Equation for vy3------------------------
                

    end
end


% Showing matrix structure
figure(5)
spy(L);

% Solving the matrix
S=L\R;

end

