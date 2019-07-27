% Function gridchange()
% This function change geometry of the grid 
% based on refinement criteria
% Function return modified nodal and cell arrays
function[celnum,nodnum,celvar,celdot,celres,celpar,celnod,nodcel,nodx,nody,recnum,celrec,rennum,nodrec]=...
     gridchange(xsize,ysize,bupper,blower,bleft,bright,reslev,nstep,S,celmax,dvxmax,dvymax,drhomax,detamax,dvxmin,dvymin,drhomin,detamin,celnum,nodnum,celvar,celdot,celres,celpar,celnod,celrho,celeta,nodcel,nodrho,nodeta,nodx,nody,recnum,celrec,rennum,nodrec)


% Array to mark cells that need refinement
celref=zeros(celmax,1);
% Array to mark cells that need coarsening
celcor=zeros(celmax,1);

% Mark cells for coarsening based on cell gradients
for ci=1:1:celnum
    if(celdot(ci,1)==0 && celres(ci)>0)  
        yncoarse=1;
        % Calculate velocity change across the cell
        dvx=abs(S(celvar(ci,1))-S(celvar(ci,4)));
        dvy=abs(S(celvar(ci,2))-S(celvar(ci,3)));
        % Check density and viscosity contrasts across the cell
        if(dvx>dvxmin || dvy>dvymin)
            yncoarse=0;
        else
            for ni1=1:1:4
                % Cell-node density/viscosity contrast
                drho=abs(nodrho(celnod(ci,ni1))-celrho(ci));
                deta=abs(log10(nodeta(celnod(ci,ni1))/celeta(ci)));
                if(drho>drhomin)
                    yncoarse=0; % no for coarsening
                end
                % Coarsening
                if(deta>detamin)
                    yncoarse=0; % no for coarsening
                end
            end
        end
        % Check gradients for coarsening
        if(yncoarse==1)
             % Check velocity gradients for 4 surrounding nodes
             for cin=1:1:4
                ni=celnod(ci,cin);
                % Calculate velocity gradient around the node
                [epsxy,dvx,dvy]=epsilonxy(S,ni,nodx,nody,nodcel,celnod,celvar,celres,bupper,blower,bleft,bright,xsize,ysize);
                % Check gradients 
                if(dvx>dvxmin || dvy>dvymin)
                    yncoarse=0;
                end
             end
            if(yncoarse==1)
                 % Check all cells around 4 surrounding nodes
                 for ni=1:1:4
                     for nc=1:1:4
                          ci1=nodcel(celnod(ci,ni),nc);
                          if(ci1>0)
                            for ni1=1:1:4
                                % Cell-node density/viscosity contrast
                                drho=abs(nodrho(celnod(ci1,ni1))-celrho(ci1));
                                deta=abs(log10(nodeta(celnod(ci1,ni1))/celeta(ci1)));
                                if(drho>drhomin)
                                    yncoarse=0; % no for coarsening
                                end
                                % Coarsening
                                if(deta>detamin)
                                    yncoarse=0; % no for coarsening
                                end
                            end
                            % Calculate velocity change across the cell
                            dvx=abs(S(celvar(ci1,1))-S(celvar(ci1,4)));
                            dvy=abs(S(celvar(ci1,2))-S(celvar(ci1,3)));
                            if(dvx>dvxmin || dvy>dvymin)
                                yncoarse=0;
                            end
                          end
                     end
                 end
                 % Mark cell for coarsening
                 if(yncoarse==1)
                    celcor(ci)=1;
                 end
            end
        end
    end
end

% Mark cells for refinment based on cell gradients
for ci=1:1:celnum
    if(celdot(ci,1)==0)  
        % Calculate velocity change across the cell
        dvx=abs(S(celvar(ci,1))-S(celvar(ci,4)));
        dvy=abs(S(celvar(ci,2))-S(celvar(ci,3)));
        % Check density and viscosity contrasts across the cell
        densyn=0;
        viscyn=0;
        for ni1=1:1:4
            % Cell-node density/viscosity contrast
            drho=abs(nodrho(celnod(ci,ni1))-celrho(ci));
            deta=abs(log10(nodeta(celnod(ci,ni1))/celeta(ci)));
            % Refinement
            if(drho>drhomax)
                densyn=1; % yes for refinement
            end
            % Refinement
            if(deta>detamax)
                viscyn=1; % yes for refinement
            end
        end
        % Check gradients for refinement
        if(densyn>0 || viscyn>0 || dvx>dvxmax || dvy>dvymax)
             % Split all cells around 4 surrounding nodes
             for ni=1:1:4
                 for nc=1:1:4
                      ci1=nodcel(celnod(ci,ni),nc);
                      if(ci1>0)
                          if(celres(ci1)<reslev)
                           % Mark current cell
                           celref(ci1)=1;
                          end
                          if(celcor(ci1)==1)
                           % Mark current cell
                           celcor(ci1)=0;
                          end
                      end
                 end
             end
        end
    end
end



% Mark cells for refinment based on nodes gradients
for ni=1:1:nodnum
    % Calculate velocity gradient around the node
    [epsxy,dvx,dvy]=epsilonxy(S,ni,nodx,nody,nodcel,celnod,celvar,celres,bupper,blower,bleft,bright,xsize,ysize);
    % Check gradients 
    if(dvx>dvxmax || dvy>dvymax)
         % Split all cells around the node
         for nc=1:1:4
              ci1=nodcel(ni,nc);
              if(ci1>0 && celres(ci1)<reslev)
               % Mark current cell
               celref(ci1)=1;
              end
         end
    end
end





% Refinement cycle
ynref=1; %refinement needed y/n
refnum=0;
while (ynref==1 || refnum<2)
    
% Mark cells for refinment based on nearest cells
% Array to mark cells that need refinement
ynref=0;
for ci=1:1:celnum
    if(celdot(ci,1)==0)  
         % Check all cells around 4 surrounding nodes
         nftotal=0; % Total number of surrounding faces
         nfhigh=0; % Number of high resolution faces
         for ni=1:1:4
             for nc=1:1:4
                  ci1=nodcel(celnod(ci,ni),nc);
                   % Count number of faces
                   if(ni~=nc)
                       nftotal=nftotal+1; % Total number of surrounding faces
                       if(ci1==0 || (ci1>0 && celres(ci1)+celref(ci1)>celres(ci)+celref(ci)))
                           nfhigh=nfhigh+1; % Number of high-resolution faces
                       end
                   end
                  if(ci1>0 && ci1~=ci && celdot(ci1,1)==0)
                       % Mark current cell
                       if(celres(ci1)+celref(ci1)>celres(ci)+celref(ci)+1 && celres(ci)<reslev)
                        celref(ci)=1;
                        ynref=1;
                       end
                       % Mark surrounding cell
                       if(celres(ci1)+celref(ci1)<celres(ci)+celref(ci)-1 && celres(ci1)<reslev)
                            % Refinement
                            if (celref(ci1)==0)
                                celref(ci1)=1;
                                ynref=1;
                            end
                       end
                       % Check one more layer of cells
                       for ni1=1:1:4
                           for nc1=1:1:4
                              ci2=nodcel(celnod(ci1,ni1),nc1);
                              if(ci2>0 && ci2~=ci && celdot(ci2,1)==0)
                                   % Mark current cell
                                   if(celres(ci2)+celref(ci2)>celres(ci)+celref(ci)+1 && celres(ci)<reslev)
                                        % Refinement
                                        if (celref(ci)==0)
                                            celref(ci)=1;
                                            ynref=1;
                                        end
                                   end
                               end
                           end
                       end
                  end
             end
         end
       % Mark fully surrounded cells
       if(nfhigh>=6)
        celref(ci)=1;
        ynref=1;
       end
    end
end

% Check cells marked for coarsening
if(refnum==0)
% Remove coarsening for cells scheduled for refinement
for ci=1:1:celnum
    if(celcor(ci)==1 && celref(ci)==1)
        celcor(ci)=0;
    end
end
% Remove coarsening for badly located cells    
for ci=1:1:celnum
    if(celcor(ci)==1)
        % Check 4 daughter cells
        cip=celpar(ci);
        yncoarse=1;
        for cn=1:1:4
            if(celcor(celdot(cip,cn))==0 || celref(celdot(cip,cn))==1 || celdot(celdot(cip,cn))~=0 || celres(celdot(cip,cn))~=celres(ci))
                yncoarse=0;
            end
        end
        % Check cells around 4 daughter cells
        if(yncoarse==1)
            for cn=1:1:4
                % Current cell index
                ci1=celdot(cip,cn);
                ci1ref=celref(ci1);
                for ni=1:1:4
                    for nic=1:1:4
                        % Nearest cell index
                        ci2=nodcel(celnod(ci1,ni),nic);
                        if(ci2>0 && celres(ci2)+celref(ci2)>celres(ci1))
                            yncoarse=0;
                        end
                    end
                end
            end
        end
        % Remove coarsening for 4 daughter cells
        if(yncoarse==0)
            for cn=1:1:4
               celcor(celdot(cip,cn))=0;
            end
        end
    end
end

figure(2)
m=[1 2 4 3 1];
for ci=1:1:celnum
    if(celdot(ci,1)==0 && celres(ci)>-1)
    for n=1:1:5
        xn(n)=nodx(celnod(ci,m(n)));
        yn(n)=nody(celnod(ci,m(n)));
    end
    if (celref(ci)==0 && celcor(ci)==0) plot(xn/1000,yn/1000,'k -'); end
    if (celref(ci)==1 && celcor(ci)==0) plot(xn/1000,yn/1000,'r -'); end
    if (celref(ci)==0 && celcor(ci)==1) plot(xn/1000,yn/1000,'g -'); end
    if (celref(ci)==1 && celcor(ci)==1) plot(xn/1000,yn/1000,'b -'); end
    hold on
    end
end
hold off
axis ij image;

end


% Refine marked cells
for rescur=0:1:reslev-1
    cimax=celnum;
    for ci=1:1:cimax
%         if(celref(ci)==1)
        if(celref(ci)==1 && celres(ci)==rescur)
             % Check contacting cells around 4 surrounding nodes
             ynrefine=1; % refinement mark
             for ni=1:1:4
                 for nc=1:1:4
                      ci1=nodcel(celnod(ci,ni),nc);
                       % Count number of faces
                       if(ni~=nc && ci1~=ci)
                           if(ci1>0 && celres(ci1)<celres(ci))
                             ynrefine=0; % refinement mark
                           end
                       end
                 end
             end
            % Coarsen parent cell of current cell
            if(ynrefine==1)
                % Split current cell
                [celnum,nodnum,celdot,celres,celpar,celnod,nodcel,nodx,nody,recnum,celrec,rennum,nodrec]=...
                 splitcelbest(ci,celnum,nodnum,celdot,celres,celpar,celnod,nodcel,nodx,nody,recnum,celrec,rennum,nodrec);  
                celref(ci)=0;
            end
        end
    end
end
% Reset refinement array
celref=zeros(celmax,1);


% Coarse marked cells
if(refnum==0)
for rescur=reslev:-1:1
    for ci=1:1:cimax
        if(celcor(ci)==1 && celres(ci)==rescur)
%         if(celcor(ci)==1)
             % Check contacting cells around 4 surrounding nodes
             yncoarsen=1; % coarsening mark
             for ni=1:1:4
                 for nc=1:1:4
                      ci1=nodcel(celnod(ci,ni),nc);
                       % Count number of faces
                       if(ni~=nc && ci1~=ci)
                           if(ci1>0 && celres(ci1)>celres(ci))
                             yncoarsen=0; % coarsening mark
                           end
                       end
                 end
             end
            % Coarsen parent cell of current cell
            if(yncoarsen==1)
                [celcor,celdot,celnod,nodcel,recnum,celrec,rennum,nodrec]=...
                 mergecelbest(celpar(ci),celdot,celnod,nodcel,celcor,recnum,celrec,rennum,nodrec); 
            end
        end
    end
end
end
% Reset coarsening array
celcor=zeros(celmax,1);


% Increase refinment cycle counter
refnum=refnum+1;


figure(3)
m=[1 2 4 3 1];
for ci=1:1:celnum
    if(celdot(ci,1)==0 && celres(ci)>-1)
    for n=1:1:5
        xn(n)=nodx(celnod(ci,m(n)));
        yn(n)=nody(celnod(ci,m(n)));
    end
    plot(xn/1000,yn/1000,'k -'); 
    hold on
    end
end
% plot([146],[24],'. r')
hold off
axis ij image;
title(['timestep=',num2str(nstep),' refinement=',num2str(reslev),'  iteration=',num2str(refnum)])
    

pause(1)
end
% End refinement cycle
end

