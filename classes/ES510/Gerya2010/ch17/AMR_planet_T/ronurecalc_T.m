% Function ronurecalc_planet()
% This function recalculate
% density and viscosity for nodes and cells
% from markers
% Function return modified nodal and cell arrays
function[celeta,celrho,celcp,celht,celtk,nodeta,nodrho,nodkt]=...
 ronurecalc_T(marknum,celnum,nodnum,celdot,celnod,celeta,celrho,celtk,celcp,celht,nodeta,nodrho,nodkt,nodx,nody,xsize,ysize,Nx,Ny,xm,ym,tm,tkm,ETAM,RHOM,CPM,KTM,HTM)

% Interpolating density from markers to nodes
% Nodal arrays
RHON=zeros(nodnum,1); % Nodal density array
ETAN=zeros(nodnum,1); % Nodal viscosity array
KTN=zeros(nodnum,1); % Nodal thermal conductivity array
WTN=zeros(nodnum,1); % Nodal weight array
% Cell arrays
RHOC=zeros(celnum,1); % Cell density array
ETAC=zeros(celnum,1); % Cell viscosity array
CPC=zeros(celnum,1); % Cell heat capacity array
HTC=zeros(celnum,1); % Cell heat production array
TKC=zeros(celnum,1); % Cell temperature array
WTC=zeros(celnum,1); % Cell weight array
% Grid step
dxb=xsize/(Nx-1);
dyb=ysize/(Ny-1);
% Marker Cycle
for m=1:1:marknum
    if(xm(m)>=0 && xm(m)<=xsize && ym(m)>=0 && ym(m)<=ysize)
        % Marker type
        mm2=tm(m);
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
        j=fix(xm(m)/dxb)+1; % Horizontal index
        if(j<1)
            j=1;
        end
        if(j>Nx-1)
            j=Nx-1;
        end
        i=fix(ym(m)/dyb)+1; % Vertical index
        if(i<1)
            i=1;
        end
        if(i>Ny-1)
            i=Ny-1;
        end
        % Defining cell index
        ci=(j-1)*(Ny-1)+i;
        % Defining distances to upper-left node
        ddxm=(xm(m)-nodx(celnod(ci,1)))/(nodx(celnod(ci,3))-nodx(celnod(ci,1))); % Horizontal
        ddym=(ym(m)-nody(celnod(ci,1)))/(nody(celnod(ci,2))-nody(celnod(ci,1)));
        % Vertical
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
        %
        % Adding density, viscosity and weigths
        % Adding to the central cell
%         wt1=(1-abs(ddxm-0.5)*2)*(1-abs(ddym-0.5)*2);
        wt1=(1-abs(ddxm-0.5))*(1-abs(ddym-0.5));
        RHOC(ci)=RHOC(ci)+RHOM(mm2)*wt1;
        ETAC(ci)=ETAC(ci)+ETAM(mm2)*wt1;
        CPC(ci)=CPC(ci)+RHOM(mm2)*CPM(mm2)*wt1;
        TKC(ci)=TKC(ci)+tkm(m)*RHOM(mm2)*CPM(mm2)*wt1;
        HTC(ci)=HTC(ci)+HTM(mm2)*wt1;
        WTC(ci)=WTC(ci)+wt1;
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
%             ddxm=(ddxm-0.5)*2;
%         else
%             ddxm=ddxm*2;
        end
        if(ddym>0.5)
            ni1=ni1+1;
            ddym=1-ddym;
%             ddym=(ddym-0.5)*2;
%         else
%             ddym=ddym*2;
        end
        % Weight for the node
        wt1=(1-ddxm)*(1-ddym);
        RHON(celnod(ci,ni1))=RHON(celnod(ci,ni1))+RHOM(mm2)*wt1;
        ETAN(celnod(ci,ni1))=ETAN(celnod(ci,ni1))+ETAM(mm2)*wt1;
        KTN(celnod(ci,ni1))=KTN(celnod(ci,ni1))+KTM(mm2)*wt1;
        WTN(celnod(ci,ni1))=WTN(celnod(ci,ni1))+wt1;
    end
end
% Recomputing density and viscosity for nodes
for ni=1:1:nodnum
    if(WTN(ni)>0)
        nodrho(ni)=RHON(ni)/WTN(ni);
        nodeta(ni)=ETAN(ni)/WTN(ni);
        nodkt(ni)=KTN(ni)/WTN(ni);
    end
end
% Recomputing density and viscosity for cells
for ci=1:1:celnum
    if(WTC(ci)>0)
        celrho(ci)=RHOC(ci)/WTC(ci);
        celeta(ci)=ETAC(ci)/WTC(ci);
        celht(ci)=HTC(ci)/WTC(ci);
        celcp(ci)=CPC(ci)/WTC(ci)/celrho(ci);
        celtk(ci)=TKC(ci)/WTC(ci)/celrho(ci)/celcp(ci);
    end
end
end
