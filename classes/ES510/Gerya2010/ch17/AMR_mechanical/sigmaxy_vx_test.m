% Function sigmaxy_vx_test()
% This function compose vx-component of SIGMAxy for specific node
% and add it to the global matrix
% Function return modified matrix
function[L,R]=sigmaxy_vx_test(kk,ni,nodx,nody,nodcel,celnod,celvar,celres,nodeta,vkf,L,R,bupper,blower,bleft,bright,g)

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


% dvx/dy
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

