% Function heatsolve()
% This function compose matrix and right parts
% for temperature equation
% Function return solutions for Poisson SP() and Stokes S()
function[tkm,celdt]=heatsolve(marknum,timesum,timestp,varnump,celnum,celdot,celnod,celvar,celres,celeta,celrho,celcp,celht,celtk,nodx,nody,nodcel,nodeta,nodkt,bupper,blower,bleft,bright,btupper,btlower,btleft,btright,xsize,ysize,S,xm,ym,tm,tkm,Nx,Ny,RHOM,CPM,KTM)
% timestp=0;

% Define global matrices for Temperature equation
L=sparse(varnump,varnump);
R=zeros(varnump,1);

% Composing equations
for ci=1:1:celnum
    if(celdot(ci,1)==0)
        %
        %       |ciA1 |ciA2 | 
        %   ----+-----+-----+----
        %  ciL1 |  ci       | ciR1
        %   ----+     TK    +---
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
        
        % Temperature equation
        %                  TKA
        %                  qy1
        %  TKL  qx1   HRO,TK,TK0,CP,HT   qx2  TKR
        %                  qy2
        %                  TKB
        % Compose matrix
        % RHO*CP*dT/dt + dqx/dx + dqy/dy = HT
        % qx=-KT*dT/dx
        % qy=-KT*dT/dy
        %
        % dqx/dx=(qx2-qx1)/dx
        % Add matrix for qx1=-dT/dx1*KTX1
        if(nodx(celnod(ci,1))>0)
            if(ciL1>0)
                [L]=qx(kp,ciL1,ciL2,ci,ciB1,nodx,nodkt,celnod,celvar,celres,-1/dx,L);
            else
                [L]=qx(kp,ciL2,ciL2,ciA1,ci,nodx,nodkt,celnod,celvar,celres,-1/dx,L);
            end
        else
            % Boundary condition
            % Insulation: dT/dx=0
            if(btleft>=0)
                % Constant temperature: dT/dx=2*(T-btleft)/dx
                KT1=(nodkt(celnod(ci,1))+nodkt(celnod(ci,2)))/2;
                L(kp,kp)=L(kp,kp)+2/dx*KT1;
                R(kp)=R(kp)+2*btleft/dx*KT1;
            end
        end
        % Add matrix for qx2=-dT/dx2*KTX2
        if(nodx(celnod(ci,3))<xsize)
            if(ciR1>0)
                [L]=qx(kp,ci,ciB2,ciR1,ciR2,nodx,nodkt,celnod,celvar,celres,1/dx,L);
            else
                [L]=qx(kp,ciA2,ci,ciR2,ciR2,nodx,nodkt,celnod,celvar,celres,1/dx,L);
            end
        else
            % Boundary condition
            % Insulation: dT/dx=0
            if(btright>=0)
                % Constant temperature: dT/dx=2*(btright-T)/dx
                KT1=(nodkt(celnod(ci,3))+nodkt(celnod(ci,4)))/2;
                L(kp,kp)=L(kp,kp)+2/dx*KT1;
                R(kp)=R(kp)+2*btright/dx*KT1;
            end
        end
        %
        % dqy/dy=(qy2-qy1)/dy
        % Add matrix for qy1=-dT/dy1*KTY1
        if(nody(celnod(ci,1))>0)
            if(ciA1>0)
                [L]=qy(kp,ciA1,ciA2,ci,ciR1,nody,nodkt,celnod,celvar,celres,-1/dy,L);
            else
                [L]=qy(kp,ciA2,ciA2,ciL1,ci,nody,nodkt,celnod,celvar,celres,-1/dy,L);
            end
        else
            % Boundary condition
            % Insulation: dT/dy=0
            if(btupper>=0)
                % Constant temperature: dT/dy=2*(T-btupper)/dy
                KT1=(nodkt(celnod(ci,1))+nodkt(celnod(ci,3)))/2;
                L(kp,kp)=L(kp,kp)+2/dy*KT1;
                R(kp)=R(kp)+2*btupper/dy*KT1;
            end
        end
        %
        % Add matrix for qy2=-dT/dy2*KTY2
        if(nody(celnod(ci,4))<ysize)
            if(ciB1>0)
                [L]=qy(kp,ci,ciR2,ciB1,ciB2,nody,nodkt,celnod,celvar,celres,1/dy,L);
            else
                [L]=qy(kp,ciL2,ci,ciB2,ciB2,nody,nodkt,celnod,celvar,celres,1/dy,L);
            end
        else
            % Boundary condition
            % Insulation: dT/dy=0
            if(btlower>=0)
                % Constant temperature: dT/dy=2*(blower-T)/dy
                KT1=(nodkt(celnod(ci,2))+nodkt(celnod(ci,4)))/2;
                L(kp,kp)=L(kp,kp)+2/dy*KT1;
                R(kp)=R(kp)+2*btlower/dx*KT1;
            end
        end
        %
        % RHO*CP*dT/dt
        % Add matrix for dT/dt=(TK-TK0)/dt
        % Left part
        % TK*CP*RHO/dt 
        L(kp,kp)=L(kp,kp)+celcp(ci)*celrho(ci)/timestp;
        % Right part
        % HT+TK0*CP*RHO/dt 
        R(kp)=R(kp)+celtk(ci)*celcp(ci)*celrho(ci)/timestp;
        %
        % Radioactive heating
        R(kp)=R(kp)+celht(ci);
        %
        % Adiabvatic heating ????
        %
        % shear heating ????
    end
end


% % Showing matrix structure
% figure(4)
% spy(L);

% Solving the matrix
ST=L\R;

% Reload solution
celdt=zeros(celnum,1);
celtn=zeros(celnum,1);
for ci=1:1:celnum
    if(celdot(ci,1)==0)
        celtn(ci)=ST(celvar(ci,6));
        celdt(ci)=celtn(ci)-celtk(ci);
    end
end

% Subgrid diffusion on markers
dsubgridt=1;
% Cell arrays
% DTC=zeros(celnum,1); % Temperature change array
% WTC=zeros(celnum,1); % Weight array
% Grid step
dxb=xsize/(Nx-1);
dyb=ysize/(Ny-1);

% Boundary conditions for temprature changes
if(btleft>=0)
    bdleft=0;
else
    bdleft=btleft;
end
if(btright>=0)
    bdright=0;
else
    bdright=btright;
end
if(btupper>=0)
    bdupper=0;
else
    bdupper=btupper;
end
if(btlower>=0)
    bdlower=0;
else
    bdlower=btlower;
end

% Thermal iterations
if(timesum==0)
    nnnmax=1;
else
    nnnmax=15;
end
nnn=0
dtmin=min(celdt)
dtmax=max(celdt)
for nnn=1:1:nnnmax

    % Cell arrays
    DTKC=zeros(celnum,1); % Cell temperature array
    WTC=zeros(celnum,1); % Cell weight array
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
            % Match nodal temperature at the first step
            if(timesum==0)
                tkm(m)=celinter(ci,celtn,xm(m),ym(m),nodx,nody,nodcel,celnod,btleft,btright,btupper,btlower,xsize,ysize);
            else
                % Temperature change
                dtkm=0;
                % Subgrid diffusion
                if (nnn==1)
                    % Interpolating old Temperature for the marker
                    tkold=celinter(ci,celtk,xm(m),ym(m),nodx,nody,nodcel,celnod,btleft,btright,btupper,btlower,xsize,ysize);
                    % Computing old temperature difference 
                    dtkm=tkold-tkm(m);
                    % Cell sizes
                    dx=(nodx(celnod(ci,3))-nodx(celnod(ci,1))); 
                    dy=(nody(celnod(ci,2))-nody(celnod(ci,1))); 
                    % Compute local thermal diffusion timescale for the marker
                    tdm=RHOM(mm2)*CPM(mm2)/KTM(mm2)/(2/dx^2+2/dy^2);
                    % Computing subgrid diffusion
                    sdif=-dsubgridt*timestp/tdm;
                    if(sdif<-30) 
                        sdif=-30;
                    end
                    dtkm=dtkm*(1-exp(sdif));
                end
                % Interpolating New temperature change for the marker
                dtkm=dtkm+celinter(ci,celdt,xm(m),ym(m),nodx,nody,nodcel,celnod,bdleft,bdright,bdupper,bdlower,xsize,ysize);
                % Computing new temperature for the marker
                tkm(m)=tkm(m)+dtkm;
                % Adding changes to the cell center
                wt1=(1-abs(ddxm-0.5))*(1-abs(ddym-0.5));
                DTKC(ci)=DTKC(ci)+dtkm*RHOM(mm2)*CPM(mm2)*wt1;
                WTC(ci)=WTC(ci)+wt1;
            end
        end
    end
    % Recomputing density and viscosity for cells
    for ci=1:1:celnum
        if(WTC(ci)>0)
            celdt(ci)=celdt(ci)-DTKC(ci)/WTC(ci)/celrho(ci)/celcp(ci);
        end
    end   
    nnn
    dtmin=min(celdt)
    dtmax=max(celdt)
end
end

