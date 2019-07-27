% Solving 2D Stokes and continuity equations
% for viscoplastic rheology based on
% primitive variable formulation
% combined with FD-based AMR and markers
% ETA*(d2vx/dx2+d2vx/dy2)-dP/dx=0
% ETA*(d2vy/dx2+d2vy/dy2)-dP/dy=-RHO*gy
% dxv/dx+dvy/dy=0

% Clearing
clear all; clf;


% Model description
Nx=31; % horizontal resolution
Ny=11; % vertical resolution
gx=0; % Horizontal gravity
gy=9.81; % Vertical gravity

% Maximal number of higher resolution levels for AMR
resmax=4;
% Number of step
stepnumber=1000;

% Boundary conditions, Model Size
xsize=30000; % model width, m
ysize=10000;  % model high, m
g.xmin=0;
g.xmax= xsize;
g.ymin=0;
g.ymax= ysize;
%Boundary conditions for Stokes, [orthogonal velocity, parallel velocity(0=no slip, 1=free slip)]
bleft=[-3e-10 1];  %Left boundary velocity
bright=[3e-10 1]; %Right boundary velocity
bupper=[2e-10 1]; %Upper boundary velocity
blower=[0 1]; %Lower boundary velocity
bpres=1e+8; % Upper-right cell pressure, Pa

dx=xsize/(Nx-1); % horizontal step, m
dy=ysize/(Ny-1); % vertical step, m

% Define markers
Mx=(Nx-1)*2^(resmax+1); % Number of markers in horizontal direction
My=(Ny-1)*2^(resmax+1); % Number of markers in vertical direction
% Marker arrays
dxm=xsize/Mx; % Marker distances in horizontal direction
dym=ysize/My; % Marker distances in vertical direction
My=My*1.5;
xm=zeros(My*Mx,1); % x-coordinate
ym=zeros(My*Mx,1); % y-coordinate
typem=zeros(My*Mx,1);
RHOM=zeros(My*Mx,1); % density
ETAM=zeros(My*Mx,1); % Viscosity
margam=zeros(My*Mx,1); % Marker strain
% Set positions and properties of markers
marknum=0; % marker counter
for m=1:1:Mx
    for n=1:1:My
        % Add marker counter
        marknum=marknum+1;
        % Set Marker coordinates
        xm(marknum)=g.xmin+dxm/2+(m-1)*dxm;%+(rand-0.5)*dxm; % horizontal coordinate
        ym(marknum)=g.ymin-(g.ymax-g.ymin)/2+dym/2+(n-1)*dym;%+(rand-0.5)*dym; % vertical coordinate
            % Distances
            xdist=abs(xm(marknum)-(g.xmin+g.xmax)/2)/ysize;
            ydist=(ym(marknum)-g.ymin)/ysize;
            % Matrix
            typem(marknum)=2;
            ETAM(marknum)=1e+22;
            RHOM(marknum)=3000;
            if(fix(fix(ydist*20)/2)*2==fix(ydist*20))
                typem(marknum)=3;
            end
            % Inclusion
            if((xdist<0.30 && ydist>0.90))
                typem(marknum)=4;
                ETAM(marknum)=1e+20;
                RHOM(marknum)=3000;
            end
            % Air
            if(ydist<0.2)
                typem(marknum)=1;
                ETAM(marknum)=1e+19;
                RHOM(marknum)=1;
            end
    end
end
% End Set positions and properties of markers

% Coordinate vectors
x=g.xmin:dx:g.xmax;
y=g.ymin:dy:g.ymax;
xp=g.xmin+dx/2:dx:g.xmax-dx/2;
yp=g.ymin+dy/2:dy:g.ymax-dy/2;

% Max amount of nodes
nodmax=Ny*Nx*50;
% Max amount of cells
celmax=(Ny-1)*(Nx-1)*50;

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
bcval=0;
timesum=0;
for nstep=1:1:stepnumber
    
%Boundary conditions for Stokes, [orthogonal velocity, parallel velocity(0=no slip, 1=free slip)]
% bcval=bcval+1;
% if(bcval>10)
%     bcval=10;
% end
% bleft=[-3e-11*bcval 1];  %Left boundary velocity
% bright=[3e-11*bcval 1]; %Right boundary velocity
% bupper=[3e-11*bcval 1]; %Upper boundary velocity
% blower=[0 1]; %Lower boundary velocity
% bpres=1e+8; % Upper-right cell pressure, Pa
%     
    
    
    % Increase resolution level
 if(fix(nstep/10)*10==nstep)
       reslev=reslev+1;
 end

        
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
   
    

    % Recalculate properties of nodes and cells from markers
[celeta,celrho,nodeta,nodrho]=...
 ronurecalc_markers(marknum,celnum,nodnum,celres,celdot,celnod,celeta,celrho,nodcel,nodeta,nodrho,nodx,nody,xsize,ysize,Nx,Ny,xm,ym,ETAM,RHOM,g);
        
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
    [S, exx, eyy, exy, sxx, syy, sxy, pr]=matrixadd_viscoplastic(varnum,celnum,nodnum,celdot,celnod,celvar,celres,celeta,celrho,nodx,nody,nodcel,nodeta,nodrho,bupper,blower,bleft,bright,bpres,g,pscale,gx,gy);

    % Remap solutions to the regular grids
    koef=2^(reslev+1); if (koef>16); koef=16; end
    prn=zeros((Ny-1)*koef,(Nx-1)*koef);
    vxn=prn; 
    vyn=prn;
    etan=prn;
    rhon=prn;
    exxn=prn;
    sxxn=prn;
    eyyn=prn;
    syyn=prn;
    exyn=prn;
    sxyn=prn;
    eiin=prn;
    siin=prn;
    dx1=xsize/(Nx-1)/koef;
    dy1=ysize/(Ny-1)/koef;
    x1=g.xmin+dx1/2:dx1:g.xmax-dx1/2;
    y1=g.ymin+dy1/2:dy1:g.ymax-dy1/2;
    prmin=1e+23;
    prmax=-1e+23;
    vxmin=1e+23;
    vxmax=-1e+23;
    vymin=1e+23;
    vymax=-1e+23;
    carea=1/(Nx-1)/(Ny-1);
    for dd=0:reslev
    for ci=1:1:celnum
        if(celdot(ci,1)==0)
        if(celres(ci)==dd)
            wt=carea*0.25^celres(ci);
            % Min/Max P
            if(prmin>pr(ci))
                prmin=pr(ci);
            end
            if(prmax<pr(ci))
                prmax=pr(ci);
            end
            % Min/Max vx
            if(vxmin>S(celvar(ci,1)))
                vxmin=S(celvar(ci,1));
            end
            if(vxmax<S(celvar(ci,1)))
                vxmax=S(celvar(ci,1));
            end
            if(vxmin>S(celvar(ci,4)))
                vxmin=S(celvar(ci,4));
            end
            if(vxmax<S(celvar(ci,4)))
                vxmax=S(celvar(ci,4));
            end
            % Min/Max vy
            if(vymin>S(celvar(ci,2)))
                vymin=S(celvar(ci,2));
            end
            if(vymax<S(celvar(ci,2)))
                vymax=S(celvar(ci,2));
            end
            if(vymin>S(celvar(ci,3)))
                vymin=S(celvar(ci,3));
            end
            if(vymax<S(celvar(ci,3)))
                vymax=S(celvar(ci,3));
            end

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
                    prn(i,j)=pr(ci);
                    etan(i,j)=celeta(ci);
                    rhon(i,j)=celrho(ci,5);
                    exxn(i,j)=exx(ci);
                    sxxn(i,j)=sxx(ci);
                    eyyn(i,j)=eyy(ci);
                    syyn(i,j)=syy(ci);
                    exyn(i,j)=(exy(celnod(ci,1))+exy(celnod(ci,2))+exy(celnod(ci,3))+exy(celnod(ci,4)))/4;
                    sxyn(i,j)=(sxy(celnod(ci,1))+sxy(celnod(ci,2))+sxy(celnod(ci,3))+sxy(celnod(ci,4)))/4;
                    eiin(i,j)=(exxn(i,j)^2+exyn(i,j)^2)^0.5;
                    siin(i,j)=(sxxn(i,j)^2+sxyn(i,j)^2)^0.5;
                    
                    % vx1/vx3
                    if(x<(nodx(celnod(ci,1))+nodx(celnod(ci,3)))/2)
                        vxn(i,j)=S(celvar(ci,1));
                    else
                        vxn(i,j)=S(celvar(ci,4));
                    end
                    % vy2/vy3
                    if(y<(nody(celnod(ci,1))+nody(celnod(ci,2)))/2)
                        vyn(i,j)=S(celvar(ci,2));
                    else
                        vyn(i,j)=S(celvar(ci,3));
                    end
                end
            end
        end
        end
    end
    end
    
%     a=1000
%     pause(5)
        
    
    
 % Compute cell gradients
 dvxavr=0; %vx
 dvxmin=1e+30;
 dvxmax=0;
 dvxnum=0;
 dvyavr=0; %vy
 dvynum=0;
 dvymin=1e+30;
 dvymax=0;
 dpravr=0; %pr
 dprmin=1e+30;
 dprmax=0;
 dprnum=0;
 drhoavr=0; %rho
 drhomin=1e+30;
 drhomax=0;
 detaavr=0; %eta
 detamin=1e+30;
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
        if(dvxmin>dvx)
            dvxmin=dvx;
        end
        dvy=abs(S(celvar(ci,2))-S(celvar(ci,3)));
        dvyavr=dvyavr+dvy;
        dvynum=dvynum+1;
        if(dvymax<dvy)
            dvymax=dvy;
        end
        if(dvymin>dvy)
            dvymin=dvy;
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
            if(dvxmin>dvx)
                dvxmin=dvx;
            end
            dvyavr=dvyavr+dvy;
            dvynum=dvynum+1;
            if(dvymax<dvy)
                dvymax=dvy;
            end
            if(dvymin>dvy)
                dvymin=dvy;
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
            if(drho<drhomin)
                drhomin=drho; 
            end
            if(deta>detamax)
                detamax=deta; 
            end
            if(deta<detamin)
                detamin=deta; 
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
                    if(dpr<dprmin)
                        dprmin=dpr; 
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
    
    
    
    
    % Visualize solution
    figure(1);clf
  
    subplot(3,3,1)
    pcolor(x1,y1,prn); colorbar; shading interp;
    title('prn'); axis ij image; %caxis([-4 4]);
       
    subplot(3,3,2)
    pcolor(x1,y1,vxn); colorbar; shading interp;
    title('vxn'); axis ij image; %caxis([-1 1]);
        
    subplot(3,3,3)
    pcolor(x1,y1,vyn); colorbar; shading interp;
    title('vyn'); axis ij image; %caxis([-1 1]);
    
    subplot(3,3,4)
    pcolor(x1,y1,log10(etan)); colorbar; shading interp;
    title('etan'); axis ij image; %caxis([-4 4]);
   
    subplot(3,3,5)
    pcolor(x1,y1,rhon); colorbar; shading interp;
    title('rhon'); axis ij image; %caxis([-4 4]);
    
%     subplot(3,3,6)
%     pcolor(x1,y1,exxn); colorbar; shading interp;
%     title('exxn'); axis ij image; %caxis([-4 4]);
    
    subplot(3,3,6)
    pcolor(x1,y1,sxxn); colorbar; shading interp;
    title('sxxn'); axis ij image; %caxis([-4 4]);
% 
%     subplot(4,4,8)
%     pcolor(x1,y1,eyyn); colorbar; shading interp;
%     title('eyyn'); axis ij image; %caxis([-4 4]);
%     
%     subplot(4,4,9)
%     pcolor(x1,y1,syyn); colorbar; shading interp;
%     title('syyn'); axis ij image; %caxis([-4 4]);
% 
%     subplot(3,3,7)
%     pcolor(x1,y1,exyn); colorbar; shading interp;
%     title('exyn'); axis ij image; %caxis([-4 4]);
    
    subplot(3,3,7)
    pcolor(x1,y1,sxyn); colorbar; shading interp;
    title('sxyn'); axis ij image; %caxis([-4 4]);


    subplot(3,3,8)
    pcolor(x1,y1,eiin); colorbar; shading interp;
    title('eiin'); axis ij image; %caxis([-4 4]);
    
    subplot(3,3,9)
    pcolor(x1,y1,siin); colorbar; shading interp;
    title('siin'); axis ij image; %caxis([-4 4]);

%     subplot(3,3,8)
%     pcolor(x1,y1,typn); colorbar; shading interp;
%     title('typn'); axis ij image; %caxis([-4 4]);
% 
%     subplot(3,3,9)
%     pcolor(x1,y1,gamn); colorbar; shading interp;
%     title('gamn'); axis ij image; %caxis([-4 4]);

 
% Define timestep
timestp=1e+20;
dvx=vxmax-vxmin;
dvy=vymax-vymin;
dxb=xsize/(Nx-1);
dyb=ysize/(Ny-1);
maxdxy=0.00000001;
if(nstep>300) 
    maxdxy=0.02;
end
if(dvx>0 && timestp>dxb*maxdxy/dvx)
    timestp=dxb*maxdxy/dvx;
end
if(dvy>0 && timestp>dyb*maxdxy/dvy)
    timestp=dyb*maxdxy/dvy;
end
timestp
nstep
nodnum
celnum


% Move markers, compute stress, strain rate, strain, pressure and viscosity
[xm,ym,margam,ETAM]=...
     movemark_plastic(S,timestp,marknum,celres,celvar,celdot,celnod,nodcel,nodx,nody,xsize,ysize,Nx,Ny,xm,ym,typem,margam,ETAM,pr,exx,sxx,exy,sxy,bupper,blower,bleft,bright,g);
timesum=timesum+timestp;


% Visualizing marker type
% Pixel grid resolution
xresol=Mx+1;
yresol=My+1;
ngrid=1;
sxstp=xsize/(xresol-1);
systp=ysize/(yresol-1);
mx=g.xmin/1000:sxstp/1000:g.xmax/1000;
my=g.ymin/1000:systp/1000:g.ymax/1000;
% Process markers
markcom=NaN*ones(yresol,xresol);
markdis=1e+20*ones(yresol,xresol);
markgii=NaN*ones(yresol,xresol);
marketa=NaN*ones(yresol,xresol);
for m = 1:1:marknum
    % Define pixel cell
    m1=fix(xm(m)/sxstp)+1;
    m2=fix(ym(m)/systp)+1;
    if (m1<1)
        m1=1;
    end
    if (m1>xresol-1)
        m1=xresol-1;
    end
    if (m2<1)
        m2=1;
    end
    if (m2>yresol-1)
        m2=yresol-1;
    end
    % Define indexes of surrounding pixels
    m10min=m1-ngrid;
    if (m10min<1)
        m10min=1;
    end
    m10max=m1+1+ngrid;
    if (m10max>xresol)
        m10max=xresol;
    end
    m20min=m2-ngrid;
    if (m20min<1)
        m20min=1;
    end
    m20max=m2+1+ngrid;
    if (m20max>yresol)
        m20max=yresol;
    end
    % Update pixels around the marker
    for m10 = m10min:1:m10max
        for m20 = m20min:1:m20max 
            % Check distance to current pixel
            dx=xm(m)-(m10-1)*sxstp;
            dy=ym(m)-(m20-1)*systp;
            dd=(dx*dx+dy*dy)^0.5;
            if(dd<markdis(m20,m10))
                markcom(m20,m10)=typem(m);
                markgii(m20,m10)=margam(m);
                marketa(m20,m10)=ETAM(m);
                markdis(m20,m10)=dd;
            end
        end
    end
end

% Draw composition
figure(4); clf
subplot(4,1,1)
pcolor(mx,my,-markcom);
% colormap(gray);
box on;
colorbar;
shading interp;
axis ij image;
title(['Composition: Step=',num2str(nstep),' Myr=',num2str(timesum*1e-6/365.25/25/3600)]); 


 % Draw log viscosity
subplot(4,1,2);
pcolor(mx,my,log10(marketa));
% colormap(gray);
box on;
colorbar;
shading interp;
axis ij image;
title('log10(viscosity)'); 

% Draw strain
subplot(4,1,3);
pcolor(mx,my,(markgii));
% colormap(gray);
box on;
colorbar;
shading interp;
axis ij image;
title('strain'); 

 % Draw grid
subplot(4,1,4);
m=[1 2 4 3 1];
for ci=1:1:celnum
    if(celdot(ci,1)==0 && celres(ci)>-1)
        for n=1:1:5
            xn(n)=nodx(celnod(ci,m(n)));
            yn(n)=nody(celnod(ci,m(n)));
        end
        plot(xn/1000,yn/1000,'k -')
        hold on
    end
end
hold off
colorbar;
axis ij image;
title('grid'); 

% Print figure
if(fix(nstep/10)*10==nstep)
    nametif    =  ['plah' num2str(nstep)  '.tif'];
    print ('-dtiff', '-r300',nametif);
    namemat    =  ['plah' num2str(nstep)  '.mat'];
    save (namemat);
end



    % Grid modification with AMR
  if(nstep>10)
    % Criteria for cell splitting
    dvxmax1=dvxavr+(dvxmax-dvxavr)*10.50; % maximal vx contrast for splitting
    dvymax1=dvyavr+(dvymax-dvyavr)*10.50; % maximal vx contrast for splitting
    dprmax1=dpravr+(dprmax-dpravr)*10.9; % maximal pr contrast for splitting
    drhomax1=drhoavr+(drhomax-drhoavr)*10.9; % maximal density contrast for splitting
    detamax1=0.5; %detaavr+(detamax-detaavr)*0.10; % maximal log10(viscosity) contrast for splitting
    % Criteria for cell merging
    dvxmin1=-1;%dvxavr-(dvxavr-dvxmin)*0.5; % minimal vx contrast for merging
    dvymin1=-1;%dvyavr-(dvyavr-dvymin)*0.5; % minimal vx contrast for merging
    dprmin1=-1;%dpravr-(dpravr-dprmin)*(0.9); % minimal vx contrast for merging
    drhomin1=-1;%drhoavr-(drhoavr-drhomin)*(0.90); % minimal density contrast for merging
    detamin1=0.2;%detaavr-(detaavr-detamin)*0.90; % minimal log10(viscosity) contrast for merging
    [celnum,nodnum,celvar,celdot,celres,celpar,celnod,nodcel,nodx,nody,recnum,celrec,rennum,nodrec]=...
         gridchange_viscoplastic(pscale,bupper,blower,bleft,bright,reslev,nstep,S,celmax,dvxmax1,dvymax1,dprmax1,drhomax1,detamax1,dvxmin1,dvymin1,dprmin1,drhomin1,detamin1,celnum,nodnum,celvar,celdot,celres,celpar,celnod,celrho,celeta,nodcel,nodrho,nodeta,nodx,nody,recnum,celrec,rennum,nodrec,g);
  end
%   
%         % Refine cells based on cell coordinate
%         for ci=1:-1:celnum
%             % Rectangle 1
%             clx=(nodx(celnod(ci,1))+nodx(celnod(ci,3)))/2;
%             cly=(nody(celnod(ci,1))+nody(celnod(ci,2)))/2;
%             if(clx >  (g.xmin+0.2*(g.xmax - g.xmin)) && ...
%                 clx < (g.xmin+0.9*(g.xmax - g.xmin)))
%             if (cly > (g.ymin+0.1*(g.ymax - g.ymin)) && ...
%                 cly < (g.ymin+0.8*(g.ymax - g.ymin)))
%              if(clx<1.5*cly)
% 
%                 % Split current cell
%                 [celnum,nodnum,celdot,celres,celpar,celnod,nodcel,nodx,nody,recnum,celrec,rennum,nodrec]=...
%                  splitcelbest(ci,celnum,nodnum,celdot,celres,celpar,celnod,nodcel,nodx,nody,recnum,celrec,rennum,nodrec);  
%                 yn=1;
%              end
%             end
%             end
% 
%             
%         end  


end
