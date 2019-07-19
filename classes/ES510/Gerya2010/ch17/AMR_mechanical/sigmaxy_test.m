% Function sigmaxy_test()
% This function compose SIGMAxy for specific node
% and add it to the global matrix
% Function return modified matrix
function[L,R]=sigmaxy_test(kk,ni,nodx,nody,nodcel,celnod,celvar,celres,nodeta,vkf,L,R,bupper,blower,bleft,bright,g)

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
% Count number of coarse cells
celmn=0;
for ci0=1:1:4
    ci=nodcel(ni,ci0);
    if(ci>0 && celres(ci)==resmin)
        celmn=celmn+1;
    end
end

% dvy/dx
%  ci1  vx1  ci3
%        |
%  vy1---+---vy2
%        |
%  ci2  vx2  ci4
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
end
% Left boundary
if(nodx(ni)==g.xmin)
    % Free slip: dvy/dx=0
    if(bleft==0)
        % No slip: dvy/dx=2*vy2/dx
        dx=abs(nodx(celnod(ciR,cnR))-nodx(ni));
        L(kk,celvar(ciR,cvR))=L(kk,celvar(ciR,cvR))+vkf*2*nodeta(ni)/dx;
    end
    if(bleft==2)
        % Inclusion test
        dx=abs(nodx(celnod(ciR,cnR))-nodx(ni));
        L(kk,celvar(ciR,cvR))=L(kk,celvar(ciR,cvR))+vkf*nodeta(ni)/dx;
        x    = nodx(ni)-dx/2;
        y    = nody(ni);
        sol  = eval_anal_Dani( x, y, g );
        R(kk)=R(kk)+vkf*nodeta(ni)*sol.vz/dx;
    end
end
% Right boundary
if(nodx(ni)==g.xmax)
    % Free slip: dvy/dx=0
    if(bright==0)
        % No slip: dvy/dx=-2*vy1/dx
        dx=abs(nodx(ni)-nodx(celnod(ciL,cnL)));
        L(kk,celvar(ciL,cvL))=L(kk,celvar(ciL,cvL))-vkf*2*nodeta(ni)/dx;
    end
    if(bright==2)
        % Inclusion test
        dx=abs(nodx(ni)-nodx(celnod(ciL,cnL)));
        L(kk,celvar(ciL,cvL))=L(kk,celvar(ciL,cvL))-vkf*nodeta(ni)/dx;
        x    = nodx(ni)+dx/2;
        y    = nody(ni);
        sol  = eval_anal_Dani( x, y, g );
        R(kk)=R(kk)-vkf*nodeta(ni)*sol.vz/dx;
    end
end
% Internal nodes
if(nodx(ni)>g.xmin && nodx(ni)<g.xmax)
    dx=abs(nodx(celnod(ciR,cnR))-nodx(celnod(ciL,cnL)))/2;
    L(kk,celvar(ciR,cvR))=L(kk,celvar(ciR,cvR))+vkf*nodeta(ni)/dx;
    L(kk,celvar(ciL,cvL))=L(kk,celvar(ciL,cvL))-vkf*nodeta(ni)/dx;
%     if(celres(ciL)==celres(ciR))
%         dx=(nodx(celnod(ciR,3))-nodx(celnod(ciL,1)))/2;
%         L(kk,celvar(ciR,cvR))=L(kk,celvar(ciR,cvR))+vkf*nodeta(ni)/dx;
%         L(kk,celvar(ciL,cvL))=L(kk,celvar(ciL,cvL))-vkf*nodeta(ni)/dx;
%     end
%     if(celres(ciL)<celres(ciR))
%         ciR1=nodcel(celnod(ciR,3),4);
%         dx=(nodx(celnod(ciR1,3))-nodx(celnod(ciL,1)))/2;
%         L(kk,celvar(ciR1,cvR))=L(kk,celvar(ciR1,cvR))+vkf*nodeta(ni)/dx*6/15;
%         L(kk,celvar(ciR,cvR))=L(kk,celvar(ciR,cvR))+vkf*nodeta(ni)/dx*10/15;
%         L(kk,celvar(ciL,cvL))=L(kk,celvar(ciL,cvL))-vkf*nodeta(ni)/dx*16/15;
%     end
%     if(celres(ciL)>celres(ciR))
%         ciL1=nodcel(celnod(ciL,1),2);
%         dx=(nodx(celnod(ciR,3))-nodx(celnod(ciL1,1)))/2;
%         L(kk,celvar(ciR,cvR))=L(kk,celvar(ciR,cvR))+vkf*nodeta(ni)/dx*16/15;
%         L(kk,celvar(ciL,cvL))=L(kk,celvar(ciL,cvL))-vkf*nodeta(ni)/dx*10/15;
%         L(kk,celvar(ciL1,cvL))=L(kk,celvar(ciL1,cvL))-vkf*nodeta(ni)/dx*6/15;
%     end
end
%


%  ci1  vx1  ci3
%        |
%  vy1---+---vy2
%        |
%  ci2  vx2  ci4
%
% dvx/dy
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
% Upper boundary
if(nody(ni)==g.ymin)
    % Free slip: dvx/dy=0
    if(bupper==0)
        % No slip: dvx/dy=2*vx2/dy
        dy=abs(nody(celnod(ciB,cnB))-nody(ni));
        L(kk,celvar(ciB,cvB))=L(kk,celvar(ciB,cvB))+vkf*2*nodeta(ni)/dy;
    end
    if(bupper==2)
        % Inclusion test
        dy=abs(nody(celnod(ciB,cnB))-nody(ni));
        L(kk,celvar(ciB,cvB))=L(kk,celvar(ciB,cvB))+vkf*nodeta(ni)/dy;
        x    = nodx(ni);
        y    = nody(ni)-dy/2;
        sol  = eval_anal_Dani( x, y, g );
        R(kk)=R(kk)+vkf*nodeta(ni)*sol.vx/dy;
    end
end
% Lower boundary
if(nody(ni)==g.ymax)
    % Free slip: dvx/dy=0
    if(blower==0)
        % No slip: dvx/dy=-2*vx1/dy
        dy=abs(nody(ni)-nody(celnod(ciA,cnA)));
        L(kk,celvar(ciA,cvA))=L(kk,celvar(ciA,cvA))-vkf*2*nodeta(ni)/dy;
    end
    if(blower==2)
        % Inclusion test
        dy=abs(nody(ni)-nody(celnod(ciA,cnA)));
        L(kk,celvar(ciA,cvA))=L(kk,celvar(ciA,cvA))-vkf*nodeta(ni)/dy;
        x    = nodx(ni);
        y    = nody(ni)+dy/2;
        sol  = eval_anal_Dani( x, y, g );
        R(kk)=R(kk)-vkf*nodeta(ni)*sol.vx/dy;
    end
end
% Internal nodes
if(nody(ni)>g.ymin && nody(ni)<g.ymax)
    dy=abs(nody(celnod(ciB,cnB))-nody(celnod(ciA,cnA)))/2;
    L(kk,celvar(ciB,cvB))=L(kk,celvar(ciB,cvB))+vkf*nodeta(ni)/dy;
    L(kk,celvar(ciA,cvA))=L(kk,celvar(ciA,cvA))-vkf*nodeta(ni)/dy;
%     if(celres(ciA)==celres(ciB))
%         dy=(nody(celnod(ciB,2))-nody(celnod(ciA,1)))/2;
%         L(kk,celvar(ciB,cvB))=L(kk,celvar(ciB,cvB))+vkf*nodeta(ni)/dy;
%         L(kk,celvar(ciA,cvA))=L(kk,celvar(ciA,cvA))-vkf*nodeta(ni)/dy;
%     end
%     if(celres(ciA)<celres(ciB))
%        ciB1=nodcel(celnod(ciB,2),4);
%        dy=(nody(celnod(ciB1,2))-nody(celnod(ciA,1)))/2;
%        L(kk,celvar(ciA,cvA))=L(kk,celvar(ciA,cvA))-vkf*nodeta(ni)/dy*16/15;
%        L(kk,celvar(ciB,cvB))=L(kk,celvar(ciB,cvB))+vkf*nodeta(ni)/dy*10/15;
%        L(kk,celvar(ciB1,cvB))=L(kk,celvar(ciB1,cvB))+vkf*nodeta(ni)/dy*6/15;
%     end
%     if(celres(ciA)>celres(ciB))
%        ciA1=nodcel(celnod(ciA,1),3);
%        dy=(nody(celnod(ciB,2))-nody(celnod(ciA1,1)))/2;
%        L(kk,celvar(ciA,cvA))=L(kk,celvar(ciA,cvA))-vkf*nodeta(ni)/dy*10/15;
%        L(kk,celvar(ciA1,cvA))=L(kk,celvar(ciA1,cvA))-vkf*nodeta(ni)/dy*6/15;
%        L(kk,celvar(ciB,cvB))=L(kk,celvar(ciB,cvB))+vkf*nodeta(ni)/dy*16/15;
%     end
end
end

