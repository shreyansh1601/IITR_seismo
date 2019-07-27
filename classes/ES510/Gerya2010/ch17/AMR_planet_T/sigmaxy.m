% Function sigmaxy()
% This function compose SIGMAxy for specific node
% and add it to the global matrix
% Function return modified matrix
function[L]=sigmaxy(kk,ni,nodx,nody,nodcel,celnod,celvar,celres,nodeta,vkf,L,bupper,blower,bleft,bright,xsize,ysize)

%  ci1  vx1  ci3
%        |
%  vy1---+---vy2
%        |
%  ci2  vx2  ci4
% 
% Indexes for 4 cells
ci1=nodcel(ni,1); % upper left
ci2=nodcel(ni,2); % lower left
ci3=nodcel(ni,3); % upper right
ci4=nodcel(ni,4); % lower right
% Maximal distance to the left
ciL=ci1; 
cvL=3;
cnL=2;
if(ciL==0 || (ci2>0 && celres(ciL)>celres(ci2)))
    ciL=ci2;
    cvL=2;
    cnL=1;
end
% Maximal distance to the right
ciR=ci3; 
cvR=3;
cnR=4;
if(ciR==0 || (ci4>0 && celres(ciR)>celres(ci4)))
    ciR=ci4;
    cvR=2;
    cnR=3;
end
% Maximal distance to above
ciA=ci1; 
cvA=4;
cnA=3;
if(ciA==0 || (ci3>0 && celres(ciA)>celres(ci3)))
    ciA=ci3;
    cvA=1;
    cnA=1;
end
% Maximal distance to below
ciB=ci2; 
cvB=4;
cnB=4;
if(ciB==0 || (ci4>0 && celres(ciB)>celres(ci4)))
    ciB=ci4;
    cvB=1;
    cnB=2;
end%
% dvy/dx
% Left boundary
if(nodx(ni)==0)
    % Free slip: dvy/dx=0
    if(bleft==0)
        % No slip: dvy/dx=2*vy2/dx
        dx=nodx(celnod(ciR,cnR))-nodx(ni);
        L(kk,celvar(ciR,cvR))=L(kk,celvar(ciR,cvR))+vkf*2*nodeta(ni)/dx;
    end
end
% Right boundary
if(nodx(ni)==xsize)
    % Free slip: dvy/dx=0
    if(bright==0)
        % No slip: dvy/dx=-2*vy1/dx
        dx=nodx(ni)-nodx(celnod(ciL,cnL));
        L(kk,celvar(ciL,cvL))=L(kk,celvar(ciL,cvL))-vkf*2*nodeta(ni)/dx;
    end
end
% Internal nodes
if(nodx(ni)>0 && nodx(ni)<xsize)
    dx=(nodx(celnod(ciR,cnR))-nodx(celnod(ciL,cnL)))/2;
    L(kk,celvar(ciR,cvR))=L(kk,celvar(ciR,cvR))+vkf*nodeta(ni)/dx;
    L(kk,celvar(ciL,cvL))=L(kk,celvar(ciL,cvL))-vkf*nodeta(ni)/dx;
end
%
%  ci1  vx1  ci3
%        |
%  vy1---+---vy2
%        |
%  ci2  vx2  ci4
%
% dvx/dy
% Upper boundary
if(nody(ni)==0)
    % Free slip: dvx/dy=0
    if(bupper==0)
        % No slip: dvx/dy=2*vx2/dy
        dy=nody(celnod(ciB,cnB))-nody(ni);
        L(kk,celvar(ciB,cvB))=L(kk,celvar(ciB,cvB))+vkf*2*nodeta(ni)/dy;
    end
end
% Lower boundary
if(nody(ni)==ysize)
    % Free slip: dvx/dy=0
    if(blower==0)
        % No slip: dvx/dy=-2*vx1/dy
        dy=nody(ni)-nody(celnod(ciA,cnA));
        L(kk,celvar(ciA,cvA))=L(kk,celvar(ciA,cvA))-vkf*2*nodeta(ni)/dy;
    end
end
% Internal nodes
if(nody(ni)>0 && nody(ni)<ysize)
        dy=(nody(celnod(ciB,cnB))-nody(celnod(ciA,cnA)))/2;
        L(kk,celvar(ciB,cvB))=L(kk,celvar(ciB,cvB))+vkf*nodeta(ni)/dy;
        L(kk,celvar(ciA,cvA))=L(kk,celvar(ciA,cvA))-vkf*nodeta(ni)/dy;
end
end

