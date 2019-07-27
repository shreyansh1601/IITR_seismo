% Function qx()
% This function compose qx=-KT*dT/dx for specific cell boundary
% and add it to the global matrix
% Function return modified matrix
function[L]=qx(kp,ciL1,ciL2,ciR1,ciR2,nodx,nodkt,celnod,celvar,celres,pkf,L)


% Upper left cell center coordinates
cxL1=(nodx(celnod(ciL1,1))+nodx(celnod(ciL1,4)))/2;
% Upper right cell center coordinates
cxR1=(nodx(celnod(ciR1,1))+nodx(celnod(ciR1,4)))/2;
% x-Distance between the cells
dx=cxR1-cxL1;


% qx=-KT*dT/dx
%
% TKL1      TKR1
% qx=KT*(TKL1-TKR1)/dx
% Same resolution for the left and right cells
if(celres(ciL1)==celres(ciR1))
    % Add  TKL1, TKR1,
    KT1=(nodkt(celnod(ciR1,1))+nodkt(celnod(ciR1,2)))/2;
    L(kp,celvar(ciL1,6))=L(kp,celvar(ciL1,6))+pkf/dx*KT1;
    L(kp,celvar(ciR1,6))=L(kp,celvar(ciR1,6))-pkf/dx*KT1;
else
    % TKL1      
    %           TKR1
    % TKL2
    % qx=KT*(TKL1/2+TKL2/2-TKR1)/dx
    % Left cells < right cell
    if(celres(ciL1)>celres(ciR1))
        % Add TKR1, TKL1, TKL2
        KT1=(nodkt(celnod(ciR1,1))+nodkt(celnod(ciR1,2)))/2;
        L(kp,celvar(ciR1,6))=L(kp,celvar(ciR1,6))-pkf/dx*KT1;
        L(kp,celvar(ciL1,6))=L(kp,celvar(ciL1,6))+pkf/dx/2*KT1;
        L(kp,celvar(ciL2,6))=L(kp,celvar(ciL2,6))+pkf/dx/2*KT1;
    else
        %           TKR1      
        % TKL1
        %           TKR2
        % qx=KT*(TKL1-TKR1/2-TKR2/2)/dx
        % Left cell > right cells
        % Add TKL1, TKR1, TKR2
        KT1=(nodkt(celnod(ciL1,3))+nodkt(celnod(ciL1,4)))/2;
        L(kp,celvar(ciL1,6))=L(kp,celvar(ciL1,6))+pkf/dx*KT1;
        L(kp,celvar(ciR1,6))=L(kp,celvar(ciR1,6))-pkf/dx/2*KT1;
        L(kp,celvar(ciR2,6))=L(kp,celvar(ciR2,6))-pkf/dx/2*KT1;
    end
end
   
end

