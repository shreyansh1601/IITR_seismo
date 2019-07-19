% Solving 2D Stokes, continuity, Poisson and temperature equations
% for variable viscosity/conductivity case based on 
% primitive variable formulation
% combined with FD-based AMR
% ETA*(d2vx/dx2+d2vx/dy2)-dP/dx=0
% ETA*(d2vy/dx2+d2vy/dy2)-dP/dy=-RHO*gy
% dxv/dx+dvy/dy=0

% Clearing
clear all; clf;

% Model description
xsize=14000000; % model width, m
ysize=14000000;  % model high, m
Nx=21; % horizontal resolution
Ny=21; % vertical resolution
dx=xsize/(Nx-1); % horizontal step, m
dy=ysize/(Ny-1); % vertical step, m
%Boundary conditions for Stokes 0=No slip, 1=Free slip
bleft=1; %Left boundary
bright=1; %Right boundary
bupper=1; %Upper boundary
blower=1; %Lower boundary
%Boundary conditions for Poisson
binner=xsize/2-dx/2; % Inner boundary for Poisson
bouter=binner; % Outer boundary for Poisson

%Boundary conditions for Temperature 0=insulation, >0=constant temperature
btleft=273; %Left boundary
btright=273; %Right boundary
btupper=273; %Upper boundary
btlower=273; %Lower boundary

% Maximal number of higher resolution levels for AMR
resmax=2;
% Elapsed time
timesum=0;
% Initial timestep
timestp=0;
% Number of step
stepnumber=100;


% Coordinate vectors
x=0:dx:xsize;
y=0:dy:ysize;
xp=dx/2:dx:xsize-dx/2;
yp=dy/2:dy:ysize-dy/2;

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
nodkt=zeros(nodmax,1); % thermal conductivity in nodes
nodrec=zeros(nodmax,1); % listing of recycled nodes
rennum=0; % Number of nodes for recycling
% Cell arrays
celnod=zeros(celmax,4); % 4 nodes surrounding each cell
celres=zeros(celmax,1); % resolution level for each cell
celpar=zeros(celmax,1); % parent cell
celdot=zeros(celmax,4); % 4 doughter cells
celeta=zeros(celmax,1); % viscosity in cells
celrho=zeros(celmax,1); % density in cells
celcp=zeros(celmax,1); % heat capacity in cells
celht=zeros(celmax,1); % heat production in cells
celtk=zeros(celmax,1); % temperature in cells
celrec=zeros(celmax,1); % listing of recycled cells
recnum=0; % Number of cells for recycling


% Defining markers
Mx=(Nx-1)*2^(resmax+1); % Number of markers in horizontal direction
My=(Ny-1)*2^(resmax+1); % Number of markers in vertical direction
dxm=xsize/Mx; % Marker distances in horizontal direction
dym=ysize/My; % Marker distances in vertical direction

% Marker coordinates, density, temperature
xm=zeros(My*Mx,1); %x-coordinate
ym=zeros(My*Mx,1); %y-coordinate
tkm=zeros(My*Mx,1); %T K
tm=zeros(My*Mx,1); %type
marvx=zeros(My*Mx,1); %vx
marvy=zeros(My*Mx,1); %vy
margx=zeros(My*Mx,1); %vx
margy=zeros(My*Mx,1); %vy
RHOM=zeros(100,1); %density
ETAM=zeros(100,1); % Viscosity
CPM=zeros(100,1); % Heat capacity
KTM=zeros(100,1); % thermal conductivity

% Weak layer
ETAM(1)=1e+23;
RHOM(1)=1;
KTM(1)=3000;
CPM(1)=3.3e+6;
HTM(1)=0;
ETAM(2)=1e+23;
RHOM(2)=1;
KTM(2)=3000;
CPM(2)=3.3e+6;
HTM(2)=0;
ETAM(3)=1e+23;
RHOM(3)=1;
KTM(3)=3000;
CPM(3)=3.3e+6;
HTM(3)=0;
% Mantle 
ETAM(4)=1e+24;
RHOM(4)=4000;
KTM(4)=3;
CPM(4)=1000;
HTM(4)=2.2e-8;
ETAM(5)=1e+24;
RHOM(5)=4000;
KTM(5)=3;
CPM(5)=1000;
HTM(5)=2.2e-8;
ETAM(6)=1e+24;
RHOM(6)=4000;
KTM(6)=3;
CPM(6)=1000;
HTM(6)=2.2e-8;
% Iron
ETAM(7)=1e+20;
RHOM(7)=5000;
KTM(7)=30;
CPM(7)=1000;
HTM(7)=1e-8;

% Set positions and properties of markers
marknum=0; % marker counter
for m=1:1:Mx
    for n=1:1:My
        % Add marker counter
        marknum=marknum+1;
        % Set Marker coordinates
        xm(marknum)=dxm/2+(m-1)*dxm+(rand-0.5)*dxm; %horizontal coordinate
        ym(marknum)=dym/2+(n-1)*dym+(rand-0.5)*dym; %vertical coordinate
        
        % Weak layer
        tm(marknum)=1;
        tkm(marknum)=273;
        nn=fix(ym(marknum)/ysize*10)+1;
        if(fix(nn/2)*2==nn)
            tm(marknum)=2;
        end
        nn=fix(xm(marknum)/xsize*10)+1;
        if(fix(nn/2)*2==nn)
            tm(marknum)=tm(marknum)+1;
        end
        
        % Planet
        cdist=((xm(marknum)-xsize/2)^2+(ym(marknum)-ysize/2)^2)^0.5;
        if(cdist<5000000)
            tm(marknum)=4;
            tkm(marknum)=273+(5000000-cdist)/5000000*1000;
            nn=fix(ym(marknum)/ysize*10)+1;
            if(fix(nn/2)*2==nn)
                tm(marknum)=5;
            end
            nn=fix(xm(marknum)/xsize*10)+1;
            if(fix(nn/2)*2==nn)
                tm(marknum)=tm(marknum)+1;
            end
        end
        % Core
        if(((xm(marknum)-0.7*xsize/2)^2+(ym(marknum)-0.7*ysize/2)^2)^0.5<1000000)
            tm(marknum)=7;
            tkm(marknum)=2000;
        end
    end
end


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
if(reslev>resmax)
    reslev=resmax;
end

% Recalculate properties of nodes and cells from markers
[celeta,celrho,celcp,celht,celtk,nodeta,nodrho,nodkt]=...
 ronurecalc_T(marknum,celnum,nodnum,celdot,celnod,celeta,celrho,celtk,celcp,celht,nodeta,nodrho,nodkt,nodx,nody,xsize,ysize,Nx,Ny,xm,ym,tm,tkm,ETAM,RHOM,CPM,KTM,HTM);



% Define indexes of unknowns for Stokes
%     vy2
% vx1 P5   vx4
%     vy3
varnum=0; % Number of unknowns for stokes and continuity
varnump=0; % Number of unknowns for poisson
celvar=zeros(celmax,6); % indexes of unknowns
for ci=1:1:(Nx-1)*(Ny-1)
    [celvar,varnum,varnump]=varindex_planet(ci,varnum,varnump,nodcel,celnod,celvar,celdot);
end



% Pressure scale
pscale=min(min(nodeta(1:nodnum)))/(dx/2+dy/2);
% Composing & solving Matrix
[S,SP,SG]=matrixadd_planet(timestp,varnum,varnump,celnum,celdot,celnod,celvar,celres,celeta,celrho,nodx,nody,nodcel,nodeta,nodrho,bupper,blower,bleft,bright,binner,bouter,xsize,ysize,pscale);



% Reload solutions
koef=2^(reslev-1);
dx=xsize/(Nx-1)/koef;
dy=ysize/(Ny-1)/koef;
pr=zeros((Nx-1)*koef+1,(Nx-1)*koef+1);
vx=pr; vy=pr; eta=pr; rho=pr; wt=pr;

i=1; j=1;
for ci=1:1:celnum
    if(celdot(ci,1)==0)
        for y=nody(celnod(ci,1)):dy:nody(celnod(ci,2))
            for x=nodx(celnod(ci,1)):dx:nodx(celnod(ci,3))
                i=fix(y/dy+0.5)+1;
                j=fix(x/dx+0.5)+1;
                pr(i,j)=pr(i,j)+S(celvar(ci,5))*pscale;
                vx(i,j)=vx(i,j)+(S(celvar(ci,1))+S(celvar(ci,4)))/2;
                vy(i,j)=vy(i,j)+(S(celvar(ci,2))+S(celvar(ci,3)))/2;
                eta(i,j)=eta(i,j)+celeta(ci);
                rho(i,j)=rho(i,j)+celrho(ci);
                wt(i,j)=wt(i,j)+1;
            end
        end
        i=i+1;
        if(i>Ny-1)
            i=1;
            j=j+1;
        end
    end
end
pr=pr./wt;
vx=vx./wt;
vy=vy./wt;
eta=eta./wt;
rho=rho./wt;



varnum
prmin=min(min(pr))
prmax=max(max(pr))
vxmin=min(min(vx))
vxmax=max(max(vx))
vymin=min(min(vy))
vymax=max(max(vy))


% Moving markers by velocity field
timestp=0;
if(nstep>resmax)
    % Criteria for stepping
    dvxmax=(vxmax-vxmin);
    dvymax=(vymax-vymin);
    timestp=1e+17;
    kfxy=1;
    if(dvxmax>0 && kfxy*dxm/dvxmax<timestp)
        timestp=kfxy*dxm/dvxmax;
    end
    if(dvymax>0 && kfxy*dym/dvymax<timestp)
        timestp=kfxy*dym/dvymax;
    end
end

% Solve temperature equation
if timestp>0
    [tkm,celdt]=heatsolve(marknum,timesum,timestp,varnump,celnum,celdot,celnod,celvar,celres,celeta,celrho,celcp,celht,celtk,nodx,nody,nodcel,nodeta,nodkt,bupper,blower,bleft,bright,btupper,btlower,btleft,btright,xsize,ysize,S,xm,ym,tm,tkm,Nx,Ny,RHOM,CPM,KTM);
end

% Move markers
[xm,ym,marvx,marvy,margx,margy]=movemark_T(timestp,marknum,celres,celvar,celdot,celnod,nodcel,nodx,nody,xsize,ysize,Nx,Ny,xm,ym,marvx,marvy,margx,margy,S,SG,bupper,blower,bleft,bright);

% Recalculate properties of nodes and cells from markers
[celeta,celrho,celcp,celht,celtk,nodeta,nodrho,nodkt]=...
 ronurecalc_T(marknum,celnum,nodnum,celdot,celnod,celeta,celrho,celtk,celcp,celht,nodeta,nodrho,nodkt,nodx,nody,xsize,ysize,Nx,Ny,xm,ym,tm,tkm,ETAM,RHOM,CPM,KTM,HTM);



% Visualization
% Visualizing marker type
% Pixel grid resolution
xresol=1*(Mx+1);
yresol=1*(My+1);
ngrid=0;
sxstp=xsize/(xresol-1);
systp=ysize/(yresol-1);
% Process markers
markcom=NaN*ones(yresol,xresol);
markvx=markcom;
markvy=markcom;
markgx=markcom;
markgy=markcom;
marktk=markcom;
markdis=1e+20*ones(yresol,xresol);
for mm1 = 1:1:marknum
if(xm(mm1)>=0 && xm(mm1)<=xsize && ym(mm1)>=0 && ym(mm1)<=ysize)
    % Define pixel cell
    m1=fix(xm(mm1)/sxstp)+1;
    m2=fix(ym(mm1)/systp)+1;
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
            mdx=xm(mm1)-(m10-1)*sxstp;
            mdy=ym(mm1)-(m20-1)*systp;
            dd=(mdx*mdx+mdy*mdy)^0.5;
            if(dd<markdis(m20,m10))
                markcom(m20,m10)=tm(mm1);
                markvx(m20,m10)=marvx(mm1);
                markvy(m20,m10)=marvy(mm1);
                markgx(m20,m10)=margx(mm1);
                markgy(m20,m10)=margy(mm1);
                marktk(m20,m10)=tkm(mm1);
                markdis(m20,m10)=dd;
            end
        end
    end
end
end


% Visualize PSI using pcolor 
figure(1)


% subplot(2,3,1)
% m=[1 2 4 3 1];
% for ci=1:1:celnum
%     if(celdot(ci,1)==0 && celres(ci)>-1)
%     for n=1:1:5
%         xn(n)=nodx(celnod(ci,m(n)));
%         yn(n)=nody(celnod(ci,m(n)));
%     end
%     if (celeta(ci)>3e+19 && celeta(ci)<3e+21) plot(xn/1000,yn/1000,'k -'); end
%     if (celeta(ci)<=3e+19) plot(xn/1000,yn/1000,'b -'); end
%     if (celeta(ci)>=3e+21) plot(xn/1000,yn/1000,'r -'); end
%     hold on
%     end
% end
% hold off
% axis ij image;
% title('grid');

subplot(2,3,1)
pcolor(markcom);  colorbar; shading interp;
axis ij image; % proportional dimensions, vertical axis down
title('composition');

subplot(2,3,2)
% pcolor(xp/1000,yp/1000,pr); shading interp; colorbar;
pcolor(log10(eta));  colorbar; shading interp;
axis ij image; % proportional dimensions, vertical axis down
title('viscosity');

subplot(2,3,3)
% pcolor(xp/1000,yp/1000,pr); shading interp; colorbar;
pcolor(rho); colorbar; shading interp;
axis ij image; % proportional dimensions, vertical axis down
title('density');

% subplot(2,3,4)
% % pcolor(xp/1000,yp/1000,pr); shading interp; colorbar;
% pcolor(pr); colorbar; shading interp;
% axis ij image; % proportional dimensions, vertical axis down
% title('pressure');

subplot(2,3,4)
pcolor(marktk); colorbar; shading interp;
axis ij image; % proportional dimensions, vertical axis down
title('temperature');

% subplot(2,3,5)
% % pcolor(x/1000,yp/1000,vx); shading interp; colorbar;
% pcolor(vx); colorbar; shading interp;
% axis ij image; % proportional dimensions, vertical axis down
% title('vx');

subplot(2,3,5)
pcolor((markgx.^2+markgy.^2).^0.5); colorbar; shading interp;
axis ij image; % proportional dimensions, vertical axis down
title('gravity');

% subplot(2,3,6)
% % pcolor(xp/1000,y/1000,vy); shading interp; colorbar;
% pcolor(vy); colorbar; shading interp;
% axis ij image; % proportional dimensions, vertical axis down
% title('vy');

subplot(2,3,6)
pcolor((markvx.^2+markvy.^2).^0.5); colorbar; shading interp;
axis ij image; % proportional dimensions, vertical axis down
title('velocity');

% figure(7)
% pcolor(markvy); colorbar; shading interp;
% axis ij image; % proportional dimensions, vertical axis down
% title('vy');
% pause(5)



nstep
celnum
recnum
nodnum
rennum
celnum
nodnum
Nx
Ny
varnum
marknum
if(nstep>resmax)
    reslev1=resmax;
else
    reslev1=nstep-1;
end
Nx1=(Nx-1)*2^reslev1+1;
Ny1=(Ny-1)*2^reslev1+1;
varnum1=(Nx1-1)*Ny1+Nx1*(Ny1-1)+(Nx1-1)*(Ny1-1);
varratio=varnum1/varnum;
reslev1
Nx1
Ny1
varnum1
varratio




% Grid modification with AMR
% Criteria for cell splitting
dvxmax=3000*(vxmax-vxmin)/max(Nx-1,Ny-1); % maximal vx contrast for splitting
dvymax=3000*(vymax-vymin)/max(Nx-1,Ny-1); % maximal vx contrast for splitting
drhomax=100; % maximal density contrast for splitting
detamax=1; % maximal log10(viscosity) contrast for splitting
% Criteria for cell merging
dvxmin=1000*(vxmax-vxmin)/max(Nx-1,Ny-1); % minimal vx contrast for merging
dvymin=1000*(vymax-vymin)/max(Nx-1,Ny-1); % minimal vx contrast for merging
drhomin=25; % minimal density contrast for merging
detamin=0.25; % minimal log10(viscosity) contrast for merging
% Proceed with grid modification
[celnum,nodnum,celvar,celdot,celres,celpar,celnod,nodcel,nodx,nody,recnum,celrec,rennum,nodrec]=...
     gridchange(xsize,ysize,bupper,blower,bleft,bright,reslev,nstep,S,celmax,dvxmax,dvymax,drhomax,detamax,dvxmin,dvymin,drhomin,detamin,celnum,nodnum,celvar,celdot,celres,celpar,celnod,celrho,celeta,nodcel,nodrho,nodeta,nodx,nody,recnum,celrec,rennum,nodrec);

 
% Update timestep
timestp
timesum=timesum+timestp
end




