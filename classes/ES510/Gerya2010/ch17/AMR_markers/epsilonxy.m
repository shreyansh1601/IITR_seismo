% Function epsilonxy()
% This function compose EPSILONxy for specific node
% Function return EPSILONxy
function[epsxy,dvx,dvy]=epsilonxy(S,ni,nodx,nody,nodcel,celnod,celvar,celres,bupper,blower,bleft,bright,g)


%  ci1  vx1  ci3
%           |
%  vy1---+---vy2
%           |
%  ci2  vx2  ci4
% 
% Indexes for 4 cells
ci1=nodcel(ni,1); % upper left
ci2=nodcel(ni,2); % lower left
ci3=nodcel(ni,3); % upper right
ci4=nodcel(ni,4); % lower right

% Count number of cells
resmin=1000;
celnm=0;
for ci0=1:1:4
    ci=nodcel(ni,ci0);
    if(ci>0)
        celnm=celnm+1;
        if(celres(ci)<resmin)
            resmin=celres(ci);
        end
    end
end  
% Count number of fine cells
celmn=0;
for ci0=1:1:4
    ci=nodcel(ni,ci0);
    if(ci>0 && celres(ci)==resmin)
        celmn=celmn+1;
    end
end
% Processing scheme
maxdist=0;
% Node away from the boundary
if(celnm==4 && celmn<4)
   % 1:3 cells
   if(celmn==1)
       maxdist=0;
   end
   % 3:1 cells
   if(celmn==3)
       maxdist=1;
   end
    % 2:2 cells
   if(celmn==2)
       % Diagonal
       if(celres(ci1)==celres(ci4) || celres(ci2)==celres(ci3)) 
           maxdist=0;
       else
       % Parallel
           maxdist=0;
       end
   end
end
% Maximal Distance
if(maxdist==1)
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
    end
% Minimal Distance
else
    % Minimal distance to the left
    ciL=ci1; 
    cvL=3;
    cnL=2;
    if(ciL==0 || (ci2>0 && celres(ciL)<celres(ci2)))
        ciL=ci2;
        cvL=2;
        cnL=1;
    end
    % Minimal distance to the right
    ciR=ci3; 
    cvR=3;
    cnR=4;
    if(ciR==0 || (ci4>0 && celres(ciR)<celres(ci4)))
        ciR=ci4;
        cvR=2;
        cnR=3;
    end
    % Minimal distance to above
    ciA=ci1; 
    cvA=4;
    cnA=3;
    if(ciA==0 || (ci3>0 && celres(ciA)<celres(ci3)))
        ciA=ci3;
        cvA=1;
        cnA=1;
    end
    % Minimal distance to below
    ciB=ci2; 
    cvB=4;
    cnB=4;
    if(ciB==0 || (ci4>0 && celres(ciB)<celres(ci4)))
        ciB=ci4;
        cvB=1;
        cnB=2;
    end
end

% dvy/dx
dvydx=0;
dvy=0;
% Left boundary
if(nodx(ni)==g.xmin && ciR>0)
    % Free slip: dvy/dx=0
    if(bleft(2)==0)
        % No slip: dvy/dx=2*vy2/dx
        dx=nodx(celnod(ciR,cnR))-nodx(ni);
        dvy=2*S(celvar(ciR,cvR));
        dvydx=dvy/dx;
    end
end
% Right boundary
if(nodx(ni)==g.xmax && ciL>0)
    % Free slip: dvy/dx=0
    if(bright(2)==0)
        % No slip: dvy/dx=-2*vy1/dx
        dx=nodx(ni)-nodx(celnod(ciL,cnL));
        dvy=-2*S(celvar(ciL,cvL));
        dvydx=dvy/dx;
    end
end
% Internal nodes
if(nodx(ni)>g.xmin && nodx(ni)<g.xmax && ciL>0 && ciR>0)
    dx=(nodx(celnod(ciR,cnR))-nodx(celnod(ciL,cnL)))/2;
    dvy=(S(celvar(ciR,cvR))-S(celvar(ciL,cvL)));
    dvydx=dvy/dx;
end
%
% dvx/dy
dvxdy=0;
dvx=0;
% Upper boundary
if(nody(ni)==g.ymin && ciB>0)
    % Free slip: dvx/dy=0
    if(bupper(2)==0)
        % No slip: dvx/dy=2*vx2/dy
        dy=nody(celnod(ciB,cnB))-nody(ni);
        dvx=2*S(celvar(ciB,cvB));
        dvxdy=dvx/dy;
    end
end
% Lower boundary
if(nody(ni)==g.ymax && ciA>0)
    % Free slip: dvx/dy=0
    if(blower(2)==0)
        % No slip: dvx/dy=-2*vx1/dy
        dy=nody(ni)-nody(celnod(ciA,cnA));
        dvx=-2*S(celvar(ciA,cvA));
        dvxdy=dvx/dy;
    end
end
% Internal nodes
if(nody(ni)>g.ymin && nody(ni)<g.ymax && ciB>0 && ciA>0)
        dy=(nody(celnod(ciB,cnB))-nody(celnod(ciA,cnA)))/2;
        dvx=(S(celvar(ciB,cvB))-S(celvar(ciA,cvA)));
        dvxdy=dvx/dy;
end
% Compute EPSILONxy
epsxy=(dvxdy+dvydx)/2;
dvx=abs(dvx);
dvy=abs(dvy);
end

