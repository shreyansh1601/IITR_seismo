% Function qy()
% This function compose qy=-KT*dT/dy for specific cell boundary
% and add it to the global matrix
% Function return modified matrix
function[L]=qy(kp,ciA1,ciA2,ciB1,ciB2,nody,nodkt,celnod,celvar,celres,pkf,L)

% Upper left cell center coordinates
cyA1=(nody(celnod(ciA1,1))+nody(celnod(ciA1,4)))/2;
% Lower left cell center coordinates
cyB1=(nody(celnod(ciB1,1))+nody(celnod(ciB1,4)))/2;
% y-Distance between the cells
dy=cyB1-cyA1;


% qy=-KT*dT/dy
%
% TKA1      
% TKB1
% qy=KT*(TKA1-TKB1)/dy
% Same resolution for the upper and lower cells
if(celres(ciA1)==celres(ciB1))
    % Add  TKA1, TKB1,
    KT1=(nodkt(celnod(ciB1,1))+nodkt(celnod(ciB1,3)))/2;
    L(kp,celvar(ciA1,6))=L(kp,celvar(ciA1,6))+pkf/dy*KT1;
    L(kp,celvar(ciB1,6))=L(kp,celvar(ciB1,6))-pkf/dy*KT1;
else
    % TKA1     TKA2  
    %      TKB1
    % qy=KT*(TKA1/2+TKA2/2-TKB1)/dy
    % Upper cells < lower cell
    if(celres(ciA1)>celres(ciB1))
        % Add TKB1, TKA1, TKA2
        KT1=(nodkt(celnod(ciB1,1))+nodkt(celnod(ciB1,3)))/2;
        L(kp,celvar(ciB1,6))=L(kp,celvar(ciB1,6))-pkf/dy*KT1;
        L(kp,celvar(ciA1,6))=L(kp,celvar(ciA1,6))+pkf/dy/2*KT1;
        L(kp,celvar(ciA2,6))=L(kp,celvar(ciA2,6))+pkf/dy/2*KT1;
    else
        %      TKA1
        % TKB1     TKB2  
        % qy=KT*(TKA1-TKB1/2-TKB2/2)/dy
        % Upper cell > lower cells
        % Add TKA1, TKB1, TKB2
        KT1=(nodkt(celnod(ciA1,2))+nodkt(celnod(ciA1,4)))/2;
        L(kp,celvar(ciA1,6))=L(kp,celvar(ciA1,6))+pkf/dy*KT1;
        L(kp,celvar(ciB1,6))=L(kp,celvar(ciB1,6))-pkf/dy/2*KT1;
        L(kp,celvar(ciB2,6))=L(kp,celvar(ciB2,6))-pkf/dy/2*KT1;
    end
end    
end

