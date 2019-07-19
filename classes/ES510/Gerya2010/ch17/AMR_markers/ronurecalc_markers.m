% Function ronurecalc_markers()
% This function recalculate
% density and viscosity for nodes and cells
% from markers
% Function return modified nodal and cell arrays
function[celeta,celrho,nodeta,nodrho]=...
 ronurecalc_markers(marknum,celnum,nodnum,celres,celdot,celnod,celeta,celrho,nodcel,nodeta,nodrho,nodx,nody,xsize,ysize,Nx,Ny,xm,ym,ETAM,RHOM,g)

% Interpolating density from markers to nodes
% Nodal arrays
RHON=zeros(nodnum,1); % Nodal density array
ETAN=zeros(nodnum,1); % Nodal viscosity array
WTN=zeros(nodnum,1); % Nodal weight array
% Cell arrays
RHOC=zeros(celnum,5); % Cell density array
ETAC=zeros(celnum,1); % Cell viscosity array
WTC=zeros(celnum,5); % Cell weight array
% Grid step
dxb=xsize/(Nx-1);
dyb=ysize/(Ny-1);
% Marker Cycle
for m=1:1:marknum
    if(xm(m)>=g.xmin && xm(m)<=g.xmax && ym(m)>=g.ymin && ym(m)<=g.ymax)
        % [ni1]--------------[ni3]
        %   |           ^      |
        %   |           | ddmy |
        %   |      ci   v      |
        %   |<-------->(m)     |
        %   | ddmx             |
        %   |                  |
        % [ni2]--------------[ni4]
        % Defining i,j indexes for the upper-left node
        % of the basic level=0 grid
        j=abs(fix((xm(m)-g.xmin)/dxb)+1); % Horizontal index
        if(j<1)
            j=1;
        end
        if(j>Nx-1)
            j=Nx-1;
        end
        i=abs(fix((ym(m)-g.ymin)/dyb)+1); % Vertical index
        if(i<1)
            i=1;
        end
        if(i>Ny-1)
            i=Ny-1;
        end
        
        % Defining cell index
        ci=(j-1)*(Ny-1)+i;
        % Defining normalized distances to upper-left node
        ddxc=nodx(celnod(ci,3))-nodx(celnod(ci,1));
        ddyc=nody(celnod(ci,2))-nody(celnod(ci,1));
        ddxm=((xm(m)-nodx(celnod(ci,1)))/ddxc); % Horizontal
        ddym=((ym(m)-nody(celnod(ci,1)))/ddyc); % Vertical
        % Search for the daugter cell
        % 4 subcells in the cell
        % 1    3
        %   ci
        % 2    4
        while(celdot(ci)>0)
            % Subsell index
            ci1=1;
            if(ddxm>0.5)
                ci1=ci1+2;
                ddxm=(ddxm-0.5)*2;
            else
                ddxm=ddxm*2;
            end
            if(ddym>0.5)
                ci1=ci1+1;
                ddym=(ddym-0.5)*2;
            else
                ddym=ddym*2;
            end
            % Subsel number
            ci=celdot(ci,ci1);
        end
        % Cell size
        ddxc=nodx(celnod(ci,3))-nodx(celnod(ci,1));
        ddyc=nody(celnod(ci,2))-nody(celnod(ci,1));
        wt2=1/(ddxc*ddyc);
%         wt2=1;
        % Save normalized distances
        ddxm0=ddxm;
        ddym0=ddym;
        %
        % Adding density, viscosity and weigths
        % Adding to the central cell
        wt1=wt2*(1-abs(ddxm-0.5))*(1-abs(ddym-0.5));
        RHOC(ci,5)=RHOC(ci,5)+RHOM(m)*wt1;
        ETAC(ci)=ETAC(ci)+ETAM(m)*wt1;
        WTC(ci,5)=WTC(ci,5)+wt1;
        % Adding to one of the surrounding nodes
        % around the cell
        % ni1    ni3
        %     ci
        % ni2    ni4
        % Define node index
        ni1=1;
        if(ddxm>0.5)
            ni1=ni1+2;
            ddxm=1-ddxm;
        end
        if(ddym>0.5)
            ni1=ni1+1;
            ddym=1-ddym;
        end
        % Weight for the node
        wt1=wt2*(1-ddxm)*(1-ddym);
        RHON(celnod(ci,ni1))=RHON(celnod(ci,ni1))+RHOM(m)*wt1;
        ETAN(celnod(ci,ni1))=ETAN(celnod(ci,ni1))+ETAM(m)*wt1;
        WTN(celnod(ci,ni1))=WTN(celnod(ci,ni1))+wt1;
        %
        % Adding density to coarser grid velocity nodes
        % Subcells 
        % 1  3
        % 2  4
        % Node 1
        if(ni1==1)
            % Subcel 2: Update cell to the left
            ciL=nodcel(celnod(ci,2),1); 
            if(ciL>0 && celres(ciL)<celres(ci))
                RHOC(ciL,4)=RHOC(ciL,4)+RHOM(m)*wt1;
                WTC(ciL,4)=WTC(ciL,4)+wt1;
            end
            % Subcel 3: Update cell to the above
            ciA=nodcel(celnod(ci,3),1); 
            if(ciA>0 && celres(ciA)<celres(ci))
                RHOC(ciA,3)=RHOC(ciA,3)+RHOM(m)*wt1;
                WTC(ciA,3)=WTC(ciA,3)+wt1;
            end
        end
        % Node 2
        if(ni1==2)
            % Subcel 1: Update cell to the left
            ciL=nodcel(celnod(ci,1),2); 
            if(ciL>0 && celres(ciL)<celres(ci))
                RHOC(ciL,4)=RHOC(ciL,4)+RHOM(m)*wt1;
                WTC(ciL,4)=WTC(ciL,4)+wt1;
            end
            % Subcel 4: Update cell to the below
            ciB=nodcel(celnod(ci,4),2); 
            if(ciB>0 && celres(ciB)<celres(ci))
                RHOC(ciB,2)=RHOC(ciB,2)+RHOM(m)*wt1;
                WTC(ciB,2)=WTC(ciB,2)+wt1;
            end
        end
        % Node 3
        if(ni1==3)
            % Subcel 4: Update cell to the right
            ciR=nodcel(celnod(ci,4),3); 
            if(ciR>0 && celres(ciR)<celres(ci))
                RHOC(ciR,1)=RHOC(ciR,1)+RHOM(m)*wt1;
                WTC(ciR,1)=WTC(ciR,1)+wt1;
            end
            % Subcel 1: Update cell to the above
            ciA=nodcel(celnod(ci,1),3); 
            if(ciA>0 && celres(ciA)<celres(ci))
                RHOC(ciA,3)=RHOC(ciA,3)+RHOM(m)*wt1;
                WTC(ciA,3)=WTC(ciA,3)+wt1;
            end
        end
        % Node 4
        if(ni1==4)
            % Subcel 3: Update cell to the right
            ciR=nodcel(celnod(ci,3),4); 
            if(ciR>0 && celres(ciR)<celres(ci))
                RHOC(ciR,1)=RHOC(ciR,1)+RHOM(m)*wt1;
                WTC(ciR,1)=WTC(ciR,1)+wt1;
            end
            % Subcel 2: Update cell to the below
            ciB=nodcel(celnod(ci,2),4); 
            if(ciB>0 && celres(ciB)<celres(ci))
                RHOC(ciB,2)=RHOC(ciB,2)+RHOM(m)*wt1;
                WTC(ciB,2)=WTC(ciB,2)+wt1;
            end
        end
        %
        % Adding density and viscosity to finer/equal grid velocity nodes
        % Restore normalized distances
        ddxm=ddxm0;
        ddym=ddym0;
        % vx1
        if(ddxm<0.5)
            wt1=wt2*(1-ddxm)*(1-abs(ddym-0.5));
            RHOC(ci,1)=RHOC(ci,1)+RHOM(m)*wt1;
            WTC(ci,1)=WTC(ci,1)+wt1;
            % Update cell to the left
            ciL=nodcel(celnod(ci,1),2); 
            if(ciL>0 && celres(ciL)==celres(ci))
                RHOC(ciL,4)=RHOC(ciL,4)+RHOM(m)*wt1;
                WTC(ciL,4)=WTC(ciL,4)+wt1;
            end
            if(ciL>0 && celres(ciL)>celres(ci))
                RHON(celnod(ciL,4))=RHON(celnod(ciL,4))+RHOM(m)*wt1;
                ETAN(celnod(ciL,4))=ETAN(celnod(ciL,4))+ETAM(m)*wt1;
                WTN(celnod(ciL,4))=WTN(celnod(ciL,4))+wt1;
            end
        % vx4
        else
            wt1=wt2*ddxm*(1-abs(ddym-0.5));
            RHOC(ci,4)=RHOC(ci,4)+RHOM(m)*wt1;
            WTC(ci,4)=WTC(ci,4)+wt1;
            % Update cell to the right
            ciR=nodcel(celnod(ci,3),4); 
            if(ciR>0 && celres(ciR)==celres(ci))
                RHOC(ciR,1)=RHOC(ciR,1)+RHOM(m)*wt1;
                WTC(ciR,1)=WTC(ciR,1)+wt1;
            end
            if(ciR>0 && celres(ciR)>celres(ci))
                RHON(celnod(ciR,2))=RHON(celnod(ciR,2))+RHOM(m)*wt1;
                ETAN(celnod(ciR,2))=ETAN(celnod(ciR,2))+ETAM(m)*wt1;
                WTN(celnod(ciR,2))=WTN(celnod(ciR,2))+wt1;
            end
        end
        % vy2
        if(ddym<0.5)
            wt1=wt2*(1-abs(ddxm-0.5))*(1-ddym);
            RHOC(ci,2)=RHOC(ci,2)+RHOM(m)*wt1;
            WTC(ci,2)=WTC(ci,2)+wt1;
            % Update cell to the above
            ciA=nodcel(celnod(ci,1),3); 
            if(ciA>0 && celres(ciA)==celres(ci))
                RHOC(ciA,3)=RHOC(ciA,3)+RHOM(m)*wt1;
                WTC(ciA,3)=WTC(ciA,3)+wt1;
            end
            if(ciA>0 && celres(ciA)>celres(ci))
                RHON(celnod(ciA,4))=RHON(celnod(ciA,4))+RHOM(m)*wt1;
                ETAN(celnod(ciA,4))=ETAN(celnod(ciA,4))+ETAM(m)*wt1;
                WTN(celnod(ciA,4))=WTN(celnod(ciA,4))+wt1;
            end
        % vy3
        else
            wt1=wt2*(1-abs(ddxm-0.5))*ddym;
            RHOC(ci,3)=RHOC(ci,3)+RHOM(m)*wt1;
            WTC(ci,3)=WTC(ci,3)+wt1;
            % Update cell to the below
            ciB=nodcel(celnod(ci,2),4); 
            if(ciB>0 && celres(ciB)==celres(ci))
                RHOC(ciB,2)=RHOC(ciB,2)+RHOM(m)*wt1;
                WTC(ciB,2)=WTC(ciB,2)+wt1;
            end
            if(ciB>0 && celres(ciB)>celres(ci))
                RHON(celnod(ciB,3))=RHON(celnod(ciB,3))+RHOM(m)*wt1;
                ETAN(celnod(ciB,3))=ETAN(celnod(ciB,3))+ETAM(m)*wt1;
                WTN(celnod(ciB,3))=WTN(celnod(ciB,3))+wt1;
            end
        end       
    end
end

% Recomputing density and viscosity for nodes
for ni=1:1:nodnum
    if(WTN(ni)>0)
        nodrho(ni)=RHON(ni)/WTN(ni);
        nodeta(ni)=ETAN(ni)/WTN(ni);
    end
end

% Recomputing density and viscosity for cells
for ci=1:1:celnum
    for n=1:1:5
        if(WTC(ci,n)>0)
            celrho(ci,n)=RHOC(ci,n)/WTC(ci,n);
        end
    end
    if(WTC(ci,5)>0)
        celeta(ci)=ETAC(ci)/WTC(ci,5);
    end
end
       
end

