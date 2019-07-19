% Function ronurecalc_markers()
% This function recalculate
% density and viscosity for nodes and cells
% from markers
% Function return modified nodal and cell arrays
function[celeta,celrho,nodeta,nodrho]=...
    ronurecalc_test2(celnum,nodnum,celdot,celnod,celeta,celrho,nodeta,nodrho,nodx,nody,xsize,ysize,Nx,Ny,g,gridtest);

% Parameters for analytic density and viscosity function
nx    = 1;
nz    = 1;
km    = 1.6*pi;
kn    = 3  *pi;
sigma = 1;
B     = 6.9;
mc    = 1e3;
mb    = 1e6;

% Recomputing density and viscosity for nodes
for ni=1:1:nodnum
        
    % Inclusion
    if (gridtest==1)
        % Matrix
        nodeta(ni)=1;
        % Inclusion
        cdist=(nodx(ni)^2+nody(ni)^2)^0.5;
        if(cdist<0.5*g.xmax)
            nodeta(ni)=mc;
        end

    % solcx
    elseif(gridtest==2)
        nodrho(ni) = cos( nx*pi*nodx(ni) )*sin( nz*pi*nody(ni) );
        nodeta(ni) = 1;
        if nodx(ni) >= 0.5*(g.xmin + g.xmax)
            nodeta(ni) = mb;
        end

    % solkz
    elseif (gridtest==3)
        nodrho(ni) = -sigma*sin(km*nody(ni))*cos(kn*nodx(ni));
        nodeta(ni) = exp(2*B*nody(ni));
            
    end
end

% Recomputing density and viscosity for cells
for ci=1:1:celnum
    
    % Inclusion
    if (gridtest==1)
        % Matrix
        celeta(ci)=1;
        % Inclusion
        cdist=(((nodx(celnod(ci,1))+nodx(celnod(ci,3)))/2)^2+((nody(celnod(ci,1))+nody(celnod(ci,2)))/2)^2)^0.5;
        if(cdist<0.5*g.xmax)
            celeta(ci)=mc;
        end

    % solcx
    elseif(gridtest==2)
        celrho(ci,5) = cos( nx*pi*(nodx(celnod(ci,1))+nodx(celnod(ci,3)))/2 )*sin( nz*pi*(nody(celnod(ci,1))+nody(celnod(ci,2)))/2 );
        celrho(ci,1) = cos( nx*pi*nodx(celnod(ci,1)))*sin( nz*pi*(nody(celnod(ci,1))+nody(celnod(ci,2)))/2 );
        celrho(ci,4) = cos( nx*pi*nodx(celnod(ci,3)))*sin( nz*pi*(nody(celnod(ci,1))+nody(celnod(ci,2)))/2 );
        celrho(ci,2) = cos( nx*pi*(nodx(celnod(ci,1))+nodx(celnod(ci,3)))/2 )*sin( nz*pi*nody(celnod(ci,1)) );
        celrho(ci,3) = cos( nx*pi*(nodx(celnod(ci,1))+nodx(celnod(ci,3)))/2 )*sin( nz*pi*nody(celnod(ci,2)) );
        celeta(ci) = 1;
        if (nodx(celnod(ci,1))+nodx(celnod(ci,3)))/2 >= 0.5*(g.xmin + g.xmax)
            celeta(ci) = mb;
        end
            
    % solkz
    elseif (gridtest==3)
            celrho(ci,5) = -sigma*sin(km*(nody(celnod(ci,1))+nody(celnod(ci,2)))/2)*cos(kn*(nodx(celnod(ci,1))+nodx(celnod(ci,3)))/2 );
            celrho(ci,1) = -sigma*sin(km*(nody(celnod(ci,1))+nody(celnod(ci,2)))/2)*cos(kn*nodx(celnod(ci,1)));
            celrho(ci,4) = -sigma*sin(km*(nody(celnod(ci,1))+nody(celnod(ci,2)))/2)*cos(kn*nodx(celnod(ci,3)));
            celrho(ci,2) = -sigma*sin(km*nody(celnod(ci,1)))*cos(kn*(nodx(celnod(ci,1))+nodx(celnod(ci,3)))/2 );
            celrho(ci,3) = -sigma*sin(km*nody(celnod(ci,2)))*cos(kn*(nodx(celnod(ci,1))+nodx(celnod(ci,3)))/2 );
            celeta(ci) = exp(2*B*(nody(celnod(ci,1))+nody(celnod(ci,2)))/2);

    end
    
end

end
