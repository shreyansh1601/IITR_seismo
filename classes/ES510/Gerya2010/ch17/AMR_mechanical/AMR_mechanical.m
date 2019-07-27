% Solving 2D Stokes, continuity, Poisson and temperature equations
% for variable viscosity/conductivity case based on
% primitive variable formulation
% combined with FD-based AMR
% ETA*(d2vx/dx2+d2vx/dy2)-dP/dx=0
% ETA*(d2vy/dx2+d2vy/dy2)-dP/dy=-RHO*gy
% dxv/dx+dvy/dy=0

% Clearing
clear all; clf;

%  Benchmarks: 1=inclusion, 2=solcx, 3=solkz
gridtest = 2;

% Model description
Nx=21; % horizontal resolution
Ny=21; % vertical resolution
% Maximal number of higher resolution levels for AMR
resmax=5;

% Criteria for cell splitting
dvxsplit=1e+9; % maximal vx contrast for splitting
dvysplit=1e+9; % maximal vx contrast for splitting
dprsplit=1e+9; % maximal pr contrast for splitting
drhosplit=1e+9; % maximal density contrast for splitting
detasplit=0.5; % maximal log10(viscosity) contrast for splitting


% Number of steps
stepnumber=200;

% Boundary conditions, Model Size
if(gridtest==1)
    % Inclusion
    xsize=2; % model width, m
    ysize=2;  % model high, m
    g.xmin=-xsize/2;
    g.xmax= xsize/2;
    g.ymin=-ysize/2;
    g.ymax= ysize/2;
    %Boundary conditions for Stokes 0=No slip, 1=Free slip, 2=inclusion
    bleft=2;  %Left boundary
    bright=2; %Right boundary
    bupper=2; %Upper boundary
    blower=2; %Lower boundary
else
    % solcx, solkz
    xsize=1; % model width, m
    ysize=1;  % model high, m
    g.xmin=0;
    g.xmax= xsize;
    g.ymin=0;
    g.ymax= ysize;
    %Boundary conditions for Stokes 0=No slip, 1=Free slip, 2=inclusion
    bleft=1;  %Left boundary
    bright=1; %Right boundary
    bupper=1; %Upper boundary
    blower=1; %Lower boundary
end
mc=1e3;
dx=xsize/((Nx-1)/2*2); % horizontal step, m
dy=ysize/((Ny-1)/2*2); % vertical step, m




% Error arrays
L1vx=zeros(stepnumber,1);
L1vxw=zeros(stepnumber,1);
L1vy=zeros(stepnumber,1);
L1vyw=zeros(stepnumber,1);
L1pr=zeros(stepnumber,1);
L1prw=zeros(stepnumber,1);
ndofcur=zeros(stepnumber,1);
L1vxc=zeros(stepnumber,1);
L1vyc=zeros(stepnumber,1);
L1prc=zeros(stepnumber,1);
L1ndof=zeros(stepnumber,1);




% Coordinate vectors
x=g.xmin:dx:g.xmax;
y=g.ymin:dy:g.ymax;
% x(11:31)=g.xmin+10*dx:dx/2:g.xmax;
% y(11:31)=g.ymin+10*dy:dy/2:g.ymax;
xp=g.xmin+dx/2:dx:g.xmax-dx/2;
yp=g.ymin+dy/2:dy:g.ymax-dy/2;

% Max amount of nodes
nodmax=Ny*Nx*100;
% Max amount of cells
celmax=(Ny-1)*(Nx-1)*100;

% Nodal arrays
nodx=zeros(nodmax,1); % x coordinates of nodes
nody=zeros(nodmax,1); % y coordinates of nodes
nodcel=zeros(nodmax,4); % 4 cells connected to each node
nodrho=zeros(nodmax,1); % density in nodes
nodeta=zeros(nodmax,1); % viscosity in nodes
nodrec=zeros(nodmax,1); % listing of recycled nodes
rennum=0; % Number of nodes for recycling
% Cell arrays
celnod=zeros(celmax,4); % 4 nodes surrounding each cell
celres=zeros(celmax,1); % resolution level for each cell
celpar=zeros(celmax,1); % parent cell
celdot=zeros(celmax,4); % 4 doughter cells
celeta=zeros(celmax,1); % viscosity in cells
celrho=zeros(celmax,5); % density in cells
celrec=zeros(celmax,1); % listing of recycled cells
recnum=0; % Number of cells for recycling


% Resolution level=0
% Initial nodes and cells
for j=1:1:Nx
    for i=1:1:Ny
        % Node index
        ni=(j-1)*Ny+i;
        % Cell index
        ci=(j-1)*(Ny-1)+i;
        % Node coordinates
        nodx(ni)=x(j); % x coordinate
        nody(ni)=y(i); % y coordinate
        % 4 Cells connected to the node
        % 1    3
        %   ni
        % 2    4
        if(i>1 && j>1)
            nodcel(ni,1)=ci-1-(Ny-1);
        end
        if(i<Ny && j>1)
            nodcel(ni,2)=ci-(Ny-1);
        end
        if(i>1 && j<Nx)
            nodcel(ni,3)=ci-1;
        end
        if(i<Ny && j<Nx)
            nodcel(ni,4)=ci;
        end
        % 4 nodes connected to the cell
        % 1    3
        %   ci
        % 2    4
        if(i<Ny && j<Nx)
            celnod(ci,1)=ni;
            celnod(ci,2)=ni+1;
            celnod(ci,3)=ni+Ny;
            celnod(ci,4)=ni+1+Ny;
        end
    end
end
% Initial number of nodes and cells
nodnum=Nx*Ny;
celnum=(Nx-1)*(Ny-1);





% Time cycle
reslev=0; % initial resolution level
for nstep=1:1:stepnumber
    % Increase resolution level
        reslev=reslev+1;
%     if(nstep==1)
%         reslev=reslev+1;
%     end
%     if(nstep>2 && L1ndof(nstep-1)==L1ndof(nstep-2))
%         nstep 
%         L1ndof(nstep-1)
%         L1ndof(nstep-2)
%         reslev=reslev+1;
%     end
    if(reslev>resmax)
        reslev=resmax;
    end
    
  
      %Visualize grid
    figure(6)
    m=[1 2 4 3 1];
    for ci=1:1:celnum
        if(celdot(ci,1)==0 && celres(ci)>-1)
            for n=1:1:5
                xn(n)=nodx(celnod(ci,m(n)));
                yn(n)=nody(celnod(ci,m(n)));
            end
            plot(xn/1000,yn/1000,'r -')
            hold on
        end
    end
    hold off
    axis image;
    
  
    
    
    % Recalculate properties of nodes and cells from markers
    [celeta,celrho,nodeta,nodrho]=...
        ronurecalc_test2(celnum,nodnum,celdot,celnod,celeta,celrho,nodeta,nodrho,nodx,nody,xsize,ysize,Nx,Ny,g,gridtest);
    
    
    
    
    % Define indexes of unknowns for Stokes
    %     vy2
    % vx1 P5   vx4
    %     vy3
    varnum=0; % Number of unknowns for stokes and continuity
    celvar=zeros(celmax,6); % indexes of unknowns
    for ci=1:1:(Nx-1)*(Ny-1)
        [celvar,varnum]=varindex_test2(ci,varnum,nodcel,celnod,celvar,celdot,celres);
    end
    L1ndof(nstep)=varnum;
    
    
   
        
        % Recalculate properties of nodes and cells
        [celeta,celrho,nodeta,nodrho]=...
            ronurecalc_test2(celnum,nodnum,celdot,celnod,celeta,celrho,nodeta,nodrho,nodx,nody,xsize,ysize,Nx,Ny,g,gridtest);
        
        % Define indexes of unknowns for Stokes
        %     vy2
        % vx1 P5   vx4
        %     vy3
        varnum=0; % Number of unknowns for stokes and continuity
        celvar=zeros(celmax,6); % indexes of unknowns
        for ci=1:1:(Nx-1)*(Ny-1)
            [celvar,varnum]=varindex_test2(ci,varnum,nodcel,celnod,celvar,celdot,celres);
        end
        
  
    
    
    % Pressure scale
    pscale=min(min(nodeta(1:nodnum)))/(dx/2+dy/2);
    % Composing & solving Matrix
   [S]=matrixadd_test4_eta_rho(varnum,celnum,celdot,celnod,celvar,celres,celeta,celrho,nodx,nody,nodcel,nodeta,nodrho,bupper,blower,bleft,bright,g,pscale,gridtest);
%    [S]=matrixadd9(varnum,celnum,celdot,celnod,celvar,celres,celeta,celrho,nodx,nody,nodcel,nodeta,nodrho,bupper,blower,bleft,bright,g,pscale,gridtest);

    koef=2^(reslev+1);% if (koef>16); koef=16; end
    dx1=xsize/(Nx-1)/koef;
    dy1=ysize/(Ny-1)/koef;
    prn=zeros((Ny-1)*koef,(Nx-1)*koef);
    pra=prn; pre=prn;prw=prn;
    vxn=prn; vxa=prn; vxe=prn;vxw=prn;
    vyn=prn; vya=prn; vye=prn;vyw=prn;
    x1=g.xmin+dx1/2:dx1:g.xmax-dx1/2;
    y1=g.ymin+dy1/2:dy1:g.ymax-dy1/2;
    % Remap solutions to the regular grids
    carea=1/(Nx-1)/(Ny-1);
    for dd=0:reslev
    for ci=1:1:celnum
        if(celdot(ci,1)==0)
        if(celres(ci)==dd)
            wt=carea*0.25^celres(ci);
            % p5
            x=(nodx(celnod(ci,1))+nodx(celnod(ci,3)))/2;
            y=(nody(celnod(ci,1))+nody(celnod(ci,2)))/2;
            if(gridtest==1) 
                sol      = eval_anal_Dani( x, y, g );
            elseif(gridtest==2) 
                sol      = eval_sol_cx( x, y );
            elseif(gridtest==3) 
                sol      = eval_sol_kz( x, y );
            end
            pr5=sol.P;
            % vx1
            x=nodx(celnod(ci,1));
            y=(nody(celnod(ci,1))+nody(celnod(ci,2)))/2;
            if(gridtest==1) 
                sol      = eval_anal_Dani( x, y, g );
            elseif(gridtest==2) 
                sol      = eval_sol_cx( x, y );
            elseif(gridtest==3) 
                sol      = eval_sol_kz( x, y );
            end
            vx1=sol.vx;
            % vx4
            x=nodx(celnod(ci,3));
            y=(nody(celnod(ci,1))+nody(celnod(ci,2)))/2;
            if(gridtest==1) 
                sol      = eval_anal_Dani( x, y, g );
            elseif(gridtest==2) 
                sol      = eval_sol_cx( x, y );
            elseif(gridtest==3) 
                sol      = eval_sol_kz( x, y );
            end
            vx4=sol.vx;
            % vy2
            x=(nodx(celnod(ci,1))+nodx(celnod(ci,3)))/2;
            y=nody(celnod(ci,1));
            if(gridtest==1) 
                sol      = eval_anal_Dani( x, y, g );
            elseif(gridtest==2) 
                sol      = eval_sol_cx( x, y );
            elseif(gridtest==3) 
                sol      = eval_sol_kz( x, y );
            end
            vy2=sol.vz;
            % vy3
            x=(nodx(celnod(ci,1))+nodx(celnod(ci,3)))/2;
            y=nody(celnod(ci,2));
            if(gridtest==1) 
                sol      = eval_anal_Dani( x, y, g );
            elseif(gridtest==2) 
                sol      = eval_sol_cx( x, y );
            elseif(gridtest==3) 
                sol      = eval_sol_kz( x, y );
            end
            vy3=sol.vz;
            % Limits
            xbeg=nodx(celnod(ci,1))+dx1/4;
            xend=nodx(celnod(ci,3))-dx1/2;
            ybeg=nody(celnod(ci,1))+dy1/4;
            yend=nody(celnod(ci,2))-dy1/2;
            for y=ybeg:dy1:yend
                for x=xbeg:dx1:xend
                    % Pressure
                    i=fix((y-g.ymin)/dy1)+1;
                    j=fix((x-g.xmin)/dx1)+1;
                    prn(i,j)=S(celvar(ci,5))*pscale;
                    pra(i,j)=pr5;
                    pre(i,j)=prn(i,j)-pra(i,j);
                    prw(i,j)=pre(i,j)*wt;
                    % vx1/vx3
                    if(x<(nodx(celnod(ci,1))+nodx(celnod(ci,3)))/2)
                        vxn(i,j)=S(celvar(ci,1));
                        vxa(i,j)=vx1;
                    else
                        vxn(i,j)=S(celvar(ci,4));
                        vxa(i,j)=vx4;
                    end
                    vxe(i,j)=vxn(i,j)-vxa(i,j);
                    vxw(i,j)=vxe(i,j)*wt;
                    % vy2/vy3
                    if(y<(nody(celnod(ci,1))+nody(celnod(ci,2)))/2)
                        vyn(i,j)=S(celvar(ci,2));
                        vya(i,j)=vy2;
                    else
                        vyn(i,j)=S(celvar(ci,3));
                        vya(i,j)=vy3;
                    end
                    vye(i,j)=vyn(i,j)-vya(i,j);
                    vyw(i,j)=vye(i,j)*wt;
                end
            end
        end
        end
    end
    end
    
%     a=1000
%     pause(5)
        
    % Compute error
    carea=1/(Nx-1)/(Ny-1);
    SE=S.*0.0;
    prerrmax=0;
    vxerrmax=0;
    vyerrmax=0;
    celltotal=0;
    prerravr=0;
    vxerravr=0;
    vyerravr=0;
    prmin=1e+23;
    prmax=-1e+23;
    vxmin=1e+23;
    vxmax=-1e+23;
    vymin=1e+23;
    vymax=-1e+23;
    
    nnn=0;
    for ci=1:1:celnum
        if(celdot(ci,1)==0)
            % Cell area=wt
            wt=carea*0.25^celres(ci);
            % Count cells
            celltotal=celltotal+1;
            % p5
            x=(nodx(celnod(ci,1))+nodx(celnod(ci,3)))/2;
            y=(nody(celnod(ci,1))+nody(celnod(ci,2)))/2;
            if(gridtest==1) 
                sol      = eval_anal_Dani( x, y, g );
            elseif(gridtest==2) 
                sol      = eval_sol_cx( x, y );
            elseif(gridtest==3) 
                sol      = eval_sol_kz( x, y );
            end
            prerr=S(celvar(ci,5))*pscale-sol.P;
            L1pr(nstep)=L1pr(nstep)+abs(prerr)*wt;
            L1prw(nstep)=L1prw(nstep)+wt;
            % Absolute error
            SE(celvar(ci,5))=abs(prerr);
            prerravr=prerravr+abs(prerr);
            % Error integral
%             SE(celvar(ci,5))=abs(prerr)*wt;
            if(SE(celvar(ci,5))>prerrmax)
                prerrmax=SE(celvar(ci,5));
            end
            % Min/Max P
            if(prmin>S(celvar(ci,5))*pscale)
                prmin=S(celvar(ci,5))*pscale;
            end
            if(prmax<S(celvar(ci,5))*pscale)
                prmax=S(celvar(ci,5))*pscale;
            end
            % vx1
            x=nodx(celnod(ci,1));
            y=(nody(celnod(ci,1))+nody(celnod(ci,2)))/2;
            if(gridtest==1) 
                sol      = eval_anal_Dani( x, y, g );
            elseif(gridtest==2) 
                sol      = eval_sol_cx( x, y );
            elseif(gridtest==3) 
                sol      = eval_sol_kz( x, y );
            end
            vxerr=S(celvar(ci,1))-sol.vx;
            L1vx(nstep)=L1vx(nstep)+abs(vxerr)*wt/2;
            L1vxw(nstep)=L1vxw(nstep)+wt/2;
            % Absolute error
            vxerravr=vxerravr+abs(vxerr)*0.5;
            SE(celvar(ci,1))=abs(vxerr);
            % Error integral
%             SE(celvar(ci,1))=abs(vxerr)*wt;
            if(SE(celvar(ci,1))>vxerrmax)
                vxerrmax=SE(celvar(ci,1));
            end
            % Min/Max vx
            if(vxmin>S(celvar(ci,1)))
                vxmin=S(celvar(ci,1));
            end
            if(vxmax<S(celvar(ci,1)))
                vxmax=S(celvar(ci,1));
            end
            % vx4
            x=nodx(celnod(ci,3));
            y=(nody(celnod(ci,1))+nody(celnod(ci,2)))/2;
            if(gridtest==1) 
                sol      = eval_anal_Dani( x, y, g );
            elseif(gridtest==2) 
                sol      = eval_sol_cx( x, y );
            elseif(gridtest==3) 
                sol      = eval_sol_kz( x, y );
            end
            vxerr=S(celvar(ci,4))-sol.vx;
            L1vx(nstep)=L1vx(nstep)+abs(vxerr)*wt/2;
            L1vxw(nstep)=L1vxw(nstep)+wt/2;
            % Absolute error
            vxerravr=vxerravr+abs(vxerr)*0.5;
            SE(celvar(ci,4))=abs(vxerr);
            % Error integral
%             SE(celvar(ci,4))=abs(vxerr)*wt;
            if(SE(celvar(ci,4))>vxerrmax)
                vxerrmax=SE(celvar(ci,4));
            end
            % Min/Max vx
            if(vxmin>S(celvar(ci,4)))
                vxmin=S(celvar(ci,4));
            end
            if(vxmax<S(celvar(ci,4)))
                vxmax=S(celvar(ci,4));
            end
            % vy2
            x=(nodx(celnod(ci,1))+nodx(celnod(ci,3)))/2;
            y=nody(celnod(ci,1));
            if(gridtest==1) 
                sol      = eval_anal_Dani( x, y, g );
            elseif(gridtest==2) 
                sol      = eval_sol_cx( x, y );
            elseif(gridtest==3) 
                sol      = eval_sol_kz( x, y );
            end
            vyerr=S(celvar(ci,2))-sol.vz;
            L1vy(nstep)=L1vy(nstep)+abs(vyerr)*wt/2;
            L1vyw(nstep)=L1vyw(nstep)+wt/2;
            % Absolute error
            vyerravr=vyerravr+abs(vyerr)*0.5;
            SE(celvar(ci,2))=abs(vyerr);
            % Error integral
%             SE(celvar(ci,2))=abs(vyerr)*wt;
            if(SE(celvar(ci,2))>vyerrmax)
                vyerrmax=SE(celvar(ci,2));
            end
            % Min/Max vy
            if(vymin>S(celvar(ci,2)))
                vymin=S(celvar(ci,2));
            end
            if(vymax<S(celvar(ci,2)))
                vymax=S(celvar(ci,2));
            end
            % vy3
            x=(nodx(celnod(ci,1))+nodx(celnod(ci,3)))/2;
            y=nody(celnod(ci,2));
            if(gridtest==1) 
                sol      = eval_anal_Dani( x, y, g );
            elseif(gridtest==2) 
                sol      = eval_sol_cx( x, y );
            elseif(gridtest==3) 
                sol      = eval_sol_kz( x, y );
            end
            vyerr=S(celvar(ci,3))-sol.vz;
            L1vy(nstep)=L1vy(nstep)+abs(vyerr)*wt/2;
            L1vyw(nstep)=L1vyw(nstep)+wt/2;
            % Absolute error
            vyerravr=vyerravr+abs(vyerr)*0.5;
            SE(celvar(ci,3))=abs(vyerr);
            % Error integral
%             SE(celvar(ci,3))=abs(vyerr)*wt;
            if(SE(celvar(ci,3))>vyerrmax);
                vyerrmax=SE(celvar(ci,3));
            end
            % Min/Max vy
            if(vymin>S(celvar(ci,3)))
                vymin=S(celvar(ci,3));
            end
            if(vymax<S(celvar(ci,3)))
                vymax=S(celvar(ci,3));
            end
        end
    end
    prerravr=prerravr/celltotal;
    vxerravr=vxerravr/celltotal;
    vyerravr=vyerravr/celltotal;
    L1pr(nstep)=L1pr(nstep)/L1prw(nstep);
    L1vx(nstep)=L1vx(nstep)/L1vxw(nstep);
    L1vy(nstep)=L1vy(nstep)/L1vyw(nstep);
    L1prc(nstep)=prerravr;
    L1vxc(nstep)=vxerravr;
    L1vyc(nstep)=vyerravr;
    
    
 % Compute cell gradients
 dvxavr=0; %vx
 dvxmax=0;
 dvxnum=0;
 dvyavr=0; %vy
 dvynum=0;
 dvymax=0;
 dpravr=0; %pr
 dprmax=0;
 dprnum=0;
 drhoavr=0; %rho
 drhomax=0;
 detaavr=0; %eta
 detamax=0;
 detanum=0;
for ci=1:1:celnum
    if(celdot(ci,1)==0)  
        % Calculate velocity change across the cell
        dvx=abs(S(celvar(ci,1))-S(celvar(ci,4)));
        dvxavr=dvxavr+dvx;
        dvxnum=dvxnum+1;
        if(dvxmax<dvx)
            dvxmax=dvx;
        end
        dvy=abs(S(celvar(ci,2))-S(celvar(ci,3)));
        dvyavr=dvyavr+dvy;
        dvynum=dvynum+1;
        if(dvymax<dvy)
            dvymax=dvy;
        end
        % Nodes gradients
        for ni1=1:1:4
            ni=celnod(ci,ni1);
            % Calculate velocity gradient around the node
            [epsxy,dvx,dvy]=epsilonxy(S,ni,nodx,nody,nodcel,celnod,celvar,celres,bupper,blower,bleft,bright,g);
            dvxavr=dvxavr+dvx;
            dvxnum=dvxnum+1;
            if(dvxmax<dvx)
                dvxmax=dvx;
            end
            dvyavr=dvyavr+dvy;
            dvynum=dvynum+1;
            if(dvymax<dvy)
                dvymax=dvy;
            end
            % Compute viscosity and density gradients
            % Cell-node density/viscosity contrast
            drho=abs(nodrho(ni)-celrho(ci,5));
            deta=abs(log10(nodeta(ni)/celeta(ci)));
            drhoavr=drhoavr+drho;
            detaavr=detaavr+deta;
            detanum=detanum+1;
            if(drho>drhomax)
                drhomax=drho; 
            end
            if(deta>detamax)
                detamax=deta; 
            end
            % Compute pressure gradients
            for ci1=1:1:4
                ci2=nodcel(ni,ci1);
                if(ci2>0 && ci1~=ni1 && ci1~=5-ni1)
                    dpr=abs(S(celvar(ci,5))-S(celvar(ci2,5)))*pscale;
                    dpravr=dpravr+dpr;
                    dprnum=dprnum+1;
                    if(dpr>dprmax)
                        dprmax=dpr; 
                    end
                end
            end
        end
    end
end
dvxavr=dvxavr/dvxnum;
dvyavr=dvyavr/dvynum;  
dpravr=dpravr/dprnum;
drhoavr=drhoavr/detanum;
detaavr=detaavr/detanum;
prmin
prmax
vxmin
vxmax
vymin
vymax
    
    
    
    
% Plot error
figure(7); clf
subplot(2,3,1)
    plot(log10(L1ndof),log10(L1pr));
    hold on
    plot(log10(L1ndof),log10(L1prc),'r');
    hold off
    title('prerr ndof'); 
subplot(2,3,2)
    plot(log10(L1ndof),log10(L1vx));
    hold on
    plot(log10(L1ndof),log10(L1vxc),'r');
    hold off
    title('vxerr ndof'); 
subplot(2,3,3)
    plot(log10(L1ndof),log10(L1vy));
    hold on
    plot(log10(L1ndof),log10(L1vyc),'r');
    hold off
    title('vyerr ndof'); 
subplot(2,3,4)
    plot(log10(L1pr));
    hold on
    plot(log10(L1prc),'r');
    hold off
    title('prerr iter'); 
subplot(2,3,5)
    plot(log10(L1vx));
    hold on
    plot(log10(L1vxc),'r');
    hold off
    title('vxerr iter'); 
subplot(2,3,6)
    plot(log10(L1vy));
    hold on
    plot(log10(L1vyc),'r');
    hold off
    title('vyerr iter'); 
    

   
    
    
    
    % Visualize PSI using pcolor
    figure(1);clf
  
    subplot(3,3,1)
    pcolor(x1,y1,prn); colorbar; shading interp;
    title('prn'); axis image; %caxis([-4 4]);
    subplot(3,3,2)
    pcolor(x1,y1,pra); colorbar; shading interp;
    title('pra'); axis image; %caxis([-4 4]);
    subplot(3,3,3)
    pcolor(x1,y1,pre); colorbar; shading interp;
    title('pre'); axis image; 
       
    subplot(3,3,4)
    pcolor(x1,y1,vxn); colorbar; shading interp;
    title('vxn'); axis image; %caxis([-1 1]);
    subplot(3,3,5)
    pcolor(x1,y1,vxa); colorbar; shading interp;
    title('vxa'); axis image; %caxis([-1 1]);
    subplot(3,3,6)
    pcolor(x1,y1,vxe); colorbar; shading interp;
    title('vxe'); axis image; 
        
    subplot(3,3,7)
    pcolor(x1,y1,vyn); colorbar; shading interp;
    title('vyn'); axis image; %caxis([-1 1]);
    subplot(3,3,8)
    pcolor(x1,y1,vya); colorbar; shading interp;
    title('vya'); axis image; %caxis([-1 1]);
    subplot(3,3,9)
    pcolor(x1,y1,vye); colorbar; shading interp;
    title('vye'); axis image; 

% %    print ('-dtiff', '-r150','-zbuffer ','grid5'); 
%     figure(8)
%     pcolor(prn); colorbar; shading interp;
%     title('prn'); axis image; %caxis([-3 3]);
%      figure(9)
%     pcolor(pre); colorbar; shading interp;
%     title('pre'); axis image; %caxis([-.4 .4]);
    
   
    
% Grid modification with AMR
  if(nstep<stepnumber)
    % Criteria for cell splitting
    dvxmax=dvxavr+(dvxmax-dvxavr)*dvxsplit; % maximal vx contrast for splitting
    dvymax=dvyavr+(dvymax-dvyavr)*dvysplit; % maximal vx contrast for splitting
    dprmax=dpravr+(dprmax-dpravr)*dprsplit; % maximal pr contrast for splitting
    drhomax=drhoavr+(drhomax-drhoavr)*drhosplit; % maximal density contrast for splitting
    detamax=detaavr+(detamax-detaavr)*detasplit; % maximal log10(viscosity) contrast for splitting
    % Criteria for cell merging
    dvxmin=-1*(vxmax-vxmin)/max(Nx-1,Ny-1); % minimal vx contrast for merging
    dvymin=-1*(vymax-vymin)/max(Nx-1,Ny-1); % minimal vx contrast for merging
    dprmin=-1*(prmax-prmin)/max(Nx-1,Ny-1); % minimal vx contrast for merging
    drhomin=-1; % minimal density contrast for merging
    detamin=-1; % minimal log10(viscosity) contrast for merging
    % Proceed with grid modification
    prerrmax=prerravr+(prerrmax-prerravr)*10.9;
    vxerrmax=vxerravr+(vxerrmax-vxerravr)*10.9;
    vyerrmax=vyerravr+(vyerrmax-vyerravr)*10.9;
    [celnum,nodnum,celvar,celdot,celres,celpar,celnod,nodcel,nodx,nody,recnum,celrec,rennum,nodrec]=...
         gridchange_error_new(pscale,SE,prerrmax,vxerrmax,vyerrmax,bupper,blower,bleft,bright,reslev,nstep,S,celmax,dvxmax,dvymax,dprmax,drhomax,detamax,dvxmin,dvymin,dprmin,drhomin,detamin,celnum,nodnum,celvar,celdot,celres,celpar,celnod,celrho,celeta,nodcel,nodrho,nodeta,nodx,nody,recnum,celrec,rennum,nodrec,g);
  end
  
        % Refine cells based on cell coordinate
        for ci=1:-1:celnum
            % Rectangle 1
            clx=(nodx(celnod(ci,1))+nodx(celnod(ci,3)))/2;
            cly=(nody(celnod(ci,1))+nody(celnod(ci,2)))/2;
            if(clx >  (g.xmin+0.2*(g.xmax - g.xmin)) && ...
                clx < (g.xmin+0.9*(g.xmax - g.xmin)))
            if (cly > (g.ymin+0.1*(g.ymax - g.ymin)) && ...
                cly < (g.ymin+0.8*(g.ymax - g.ymin)))
             if(clx<1.5*cly)

                % Split current cell
                [celnum,nodnum,celdot,celres,celpar,celnod,nodcel,nodx,nody,recnum,celrec,rennum,nodrec]=...
                 splitcelbest(ci,celnum,nodnum,celdot,celres,celpar,celnod,nodcel,nodx,nody,recnum,celrec,rennum,nodrec);  
                yn=1;
             end
            end
            end


%             % Rectangle 2
%             if(yn==0)
%             if (nody(celnod(ci,1)) >= (g.ymin+0.45*(g.ymax - g.ymin)) && ...
%                 nody(celnod(ci,2)) <= (g.ymin+0.5*(g.ymax - g.ymin)))
%             if(nodx(celnod(ci,1)) >= (g.xmin+0.63*(g.xmax - g.xmin)) && ...
%                 nodx(celnod(ci,3)) <= (g.xmin+0.7*(g.xmax - g.xmin)))
%                 % Split current cell
%                 [celnum,nodnum,celdot,celres,celpar,celnod,nodcel,nodx,nody,recnum,celrec,rennum,nodrec]=...
%                  splitcelbest(ci,celnum,nodnum,celdot,celres,celpar,celnod,nodcel,nodx,nody,recnum,celrec,rennum,nodrec);  
%             end
%             end
%             
%             end
            
            
        end  
    
end
