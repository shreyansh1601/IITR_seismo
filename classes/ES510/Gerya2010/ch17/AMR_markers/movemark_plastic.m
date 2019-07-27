% Function movemark_T()
% This function moves markers by velocity field
% Function return modified marker arrays
function[xm,ym,margam,ETAM]=...
     movemark_plastic(S,timestp,marknum,celres,celvar,celdot,celnod,nodcel,nodx,nody,xsize,ysize,Nx,Ny,xm,ym,typem,margam,ETAM,pr,exx,sxx,exy,sxy,bupper,blower,bleft,bright,g)

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
        ddxmvx=ddxm;
        ddymvx=ddym;
        ddxmvy=ddxm;
        ddymvy=ddym;
        
        % Interpolate EPSILONii, SIGMAii
        % EPSILONxx, SIGMA'xx, Pressure
        [exx1]=celinter(ci,exx,xm(m),ym(m),nodx,nody,nodcel,celnod,g);
        [sxx1]=celinter(ci,sxx,xm(m),ym(m),nodx,nody,nodcel,celnod,g);
        [pr1]=celinter(ci,pr,xm(m),ym(m),nodx,nody,nodcel,celnod,g);
       % EPSILONxy, SIGMAxy
        exy1=(1-ddxm)*(1-ddym)*exy(celnod(ci,1))+...
            (1-ddxm)*(ddym)*exy(celnod(ci,2))+...
            (ddxm)*(1-ddym)*exy(celnod(ci,3))+...
            (ddxm)*(ddym)*exy(celnod(ci,4));
        sxy1=(1-ddxm)*(1-ddym)*sxy(celnod(ci,1))+...
            (1-ddxm)*(ddym)*sxy(celnod(ci,2))+...
            (ddxm)*(1-ddym)*sxy(celnod(ci,3))+...
            (ddxm)*(ddym)*sxy(celnod(ci,4));
        % Save EPSILONii, SIGMAii for the marker
        eii1=(exx1^2+exy1^2)^0.5;
        sii1=(sxx1^2+sxy1^2)^0.5;
        % Update strain for the marker
        margam(m)=margam(m)+eii1*timestp;
       % Compute viscosity
        % Air
        if(typem(m)==1)
            ETAM(m)=1e+19;
        end
        % Weak inclusion
        if(typem(m)==4)
            ETAM(m)=1e+20;
        end
        % Matrix
        if(typem(m)==2 || typem(m)==3)
            if(eii1>0)
                eta1=1e+22;
                % Strain weakening parameters
                strain=margam(m);
                co1=3e+7;
                co2=2e+7;
                kf1=0.6;
                kf2=0.5;
                gm1=0;
                gm2=4;
                % Strain weakening
                co=co1;
                kf=kf1;
                if(strain>gm1 && strain<gm2)
                    co=co1+(co2-co1)*(strain-gm1)/(gm2-gm1);
                    kf=kf1+(kf2-kf1)*(strain-gm1)/(gm2-gm1);
                else
                    if(strain>gm2)
                        co=co2;
                        kf=kf2;
                    end
                end
                 % Brittle/Plastic strength
                syield=co+pr1*kf;
                % Stress-based brittle/plastic strain rate
                eii1=eii1*(sii1/2/ETAM(m)/eii1)^0.93;
%                 eii1=eii1*0.5+sii1/2/ETAM(m)*0.5;
                 % Brittle/Plastic viscosity
                if(2*eta1*eii1>syield)
                    ETAM(m)=syield/2/eii1;
                    if(ETAM(m)<1e+19)
                        ETAM(m)=1e+19;
                    end
                else
                    ETAM(m)=1e+22; 
                end
            else
               ETAM(m)=1e+22; 
            end
        end
        

        
        % vx velocity interpolation
        % Interpolate down
        if(ddym>0.5)
            ddymvx=ddym-0.5;
            % Lower boundary
            if(nody(celnod(ci,2))==ysize)
                % No slip
                if(blower==0)
                    vx1=S(celvar(ci,1));
                    vx2=-S(celvar(ci,1));
                    vx3=S(celvar(ci,4));
                    vx4=-S(celvar(ci,4));
                else
                    % Free slip
                    vx1=S(celvar(ci,1));
                    vx2=S(celvar(ci,1));
                    vx3=S(celvar(ci,4));
                    vx4=S(celvar(ci,4));
                end
            else
                % Index of the lower left cell
                ci1=nodcel(celnod(ci,2),4);
                % Index of the lower right cell
                ci2=nodcel(celnod(ci,4),2);
                % Lower left cell exist
                if(ci1>0) 
                    % Same size for upper and lower cell
                    if(celres(ci)==celres(ci1))
                        vx1=S(celvar(ci,1));
                        vx2=S(celvar(ci1,1));
                        vx3=S(celvar(ci,4));
                        vx4=S(celvar(ci1,4));
                    else
                        % Recompute normalized vertical distance
                        cy=(nody(celnod(ci,2))+nody(celnod(ci,1)))/2;
                        cy1=(nody(celnod(ci1,2))+nody(celnod(ci1,1)))/2;
                        ddymvx=(ym(m)-cy)/(cy1-cy); 
                        % Upper cell > lower cell
                        if(celres(ci)<celres(ci1))
                            % Left subcell
                            if(ddxm<0.5)
                                vx1=S(celvar(ci,1));
                                vx2=S(celvar(ci1,1));
                                vx3=(S(celvar(ci,1))+S(celvar(ci,4)))/2;
                                vx4=S(celvar(ci1,4));
                                ddxmvx=ddxm*2;
                            else
                                % Right subcell
                                vx1=(S(celvar(ci,1))+S(celvar(ci,4)))/2;
                                vx2=S(celvar(ci2,1));
                                vx3=S(celvar(ci,4));
                                vx4=S(celvar(ci2,4));
                                ddxmvx=(ddxm-0.5)*2;
                            end
                        else
                            % Upper cell < lower cell
                            % Left subcell
                            vx1=S(celvar(ci,1));
                            vx2=S(celvar(ci1,1));
                            vx3=S(celvar(ci,4));
                            vx4=(S(celvar(ci1,1))+S(celvar(ci1,4)))/2;
                        end
                    end
                else
                    % Lower left cell does not exist =>
                    % Recompute normalized vertical distance
                    cy=(nody(celnod(ci,2))+nody(celnod(ci,1)))/2;
                    cy2=(nody(celnod(ci2,2))+nody(celnod(ci2,1)))/2;
                    ddymvx=(ym(m)-cy)/(cy2-cy); 
                    % => Lower right cell > upper cell
                    vx1=S(celvar(ci,1));
                    vx2=(S(celvar(ci2,1))+S(celvar(ci2,4)))/2;
                    vx3=S(celvar(ci,4));
                    vx4=S(celvar(ci2,4));
                end
            end
        else
            % Interpolate Up
            ddymvx=ddym+0.5;
            % Upper boundary
            if(nody(celnod(ci,1))==0)
                % No slip
                if(bupper==0)
                    vx1=-S(celvar(ci,1));
                    vx2=S(celvar(ci,1));
                    vx3=-S(celvar(ci,4));
                    vx4=S(celvar(ci,4));
                else
                    % Free slip
                    vx1=S(celvar(ci,1));
                    vx2=S(celvar(ci,1));
                    vx3=S(celvar(ci,4));
                    vx4=S(celvar(ci,4));
               end
            else
                % Index of the Upper left cell
                ci1=nodcel(celnod(ci,1),3);
                % Index of the Upper right cell
                ci2=nodcel(celnod(ci,3),1);
                % Upper left cell exist
                if(ci1>0) 
                    % Same size for upper and lower cell
                    if(celres(ci)==celres(ci1))
                        vx1=S(celvar(ci1,1));
                        vx2=S(celvar(ci,1));
                        vx3=S(celvar(ci1,4));
                        vx4=S(celvar(ci,4));
                    else
                         % Recompute normalized vertical distance
                        cy=(nody(celnod(ci,2))+nody(celnod(ci,1)))/2;
                        cy1=(nody(celnod(ci1,2))+nody(celnod(ci1,1)))/2;
                        ddymvx=(ym(m)-cy1)/(cy-cy1); 
                        % Upper cell < lower cell
                        if(celres(ci)<celres(ci1))
                            % Left subcell
                            if(ddxm<0.5)
                                vx1=S(celvar(ci1,1));
                                vx2=S(celvar(ci,1));
                                vx3=S(celvar(ci1,4));
                                vx4=(S(celvar(ci,1))+S(celvar(ci,4)))/2;
                                ddxmvx=ddxm*2;
                            else
                                % Right subcell
                                vx1=S(celvar(ci2,1));
                                vx2=(S(celvar(ci,1))+S(celvar(ci,4)))/2;
                                vx3=S(celvar(ci2,4));
                                vx4=S(celvar(ci,4));
                                ddxmvx=(ddxm-0.5)*2;
                            end
                        else
                            % Upper cell < lower cell
                            % Left subcell
                            vx1=S(celvar(ci1,1));
                            vx2=S(celvar(ci,1));
                            vx3=(S(celvar(ci1,1))+S(celvar(ci1,4)))/2;
                            vx4=S(celvar(ci,4));
                        end
                    end
                else
                    % Upper left cell does not exist =>
                    % Recompute normalized vertical distance
                    cy=(nody(celnod(ci,2))+nody(celnod(ci,1)))/2;
                    cy2=(nody(celnod(ci2,2))+nody(celnod(ci2,1)))/2;
                    ddymvx=(ym(m)-cy2)/(cy-cy2); 
                    % => Upper right cell > lower cell
                    vx1=(S(celvar(ci2,1))+S(celvar(ci2,4)))/2;
                    vx2=S(celvar(ci,1));
                    vx3=S(celvar(ci2,4));
                    vx4=S(celvar(ci,4));
                end
            end
        end
        % Compute vx velocity for the marker
        vxm=(1-ddxmvx)*(1-ddymvx)*vx1+...
            (1-ddxmvx)*(ddymvx)*vx2+...
            (ddxmvx)*(1-ddymvx)*vx3+...
            (ddxmvx)*(ddymvx)*vx4;

        % vy velocity interpolation
        % Interpolate right
        if(ddxm>0.5)
            ddxmvy=ddxm-0.5;
            % Right boundary
            if(nodx(celnod(ci,3))==xsize)
                % No slip
                if(bright==0)
                    vy1=S(celvar(ci,2));
                    vy2=S(celvar(ci,3));
                    vy3=-S(celvar(ci,2));
                    vy4=-S(celvar(ci,3));
                else
                    % Free slip
                    vy1=S(celvar(ci,2));
                    vy2=S(celvar(ci,3));
                    vy3=S(celvar(ci,2));
                    vy4=S(celvar(ci,3));
                end
            else
                % Index of the upper right cell
                ci1=nodcel(celnod(ci,3),4);
                % Index of the lower right cell
                ci2=nodcel(celnod(ci,4),3);
                % Lower left cell exist
                if(ci1>0) 
                    % Same size for left and right cell
                    if(celres(ci)==celres(ci1))
                        vy1=S(celvar(ci,2));
                        vy2=S(celvar(ci,3));
                        vy3=S(celvar(ci1,2));
                        vy4=S(celvar(ci1,3));
                    else
                         % Recompute normalized vertical distance
                        cx=(nodx(celnod(ci,3))+nodx(celnod(ci,1)))/2;
                        cx1=(nodx(celnod(ci1,3))+nodx(celnod(ci1,1)))/2;
                        ddxmvy=(xm(m)-cx)/(cx1-cx); 
                        % Left cell > right cell
                        if(celres(ci)<celres(ci1))
                            % Upper subcell
                            if(ddym<0.5)
                                vy1=S(celvar(ci,2));
                                vy2=(S(celvar(ci,2))+S(celvar(ci,3)))/2;
                                vy3=S(celvar(ci1,2));
                                vy4=S(celvar(ci1,3));
                                ddymvy=ddym*2;
                            else
                                % Lower subcell
                                vy1=(S(celvar(ci,2))+S(celvar(ci,3)))/2;
                                vy2=S(celvar(ci,3));
                                vy3=S(celvar(ci2,2));
                                vy4=S(celvar(ci2,3));
                                ddymvy=(ddym-0.5)*2;
                            end
                        else
                            % Left cell < Right cell
                            % Upper subcell
                            vy1=S(celvar(ci,2));
                            vy2=S(celvar(ci,3));
                            vy3=S(celvar(ci1,2));
                            vy4=(S(celvar(ci1,2))+S(celvar(ci1,3)))/2;
                        end
                    end
                else
                    % Upper right cell does not exist =>
                    % Recompute normalized vertical distance
                    cx=(nodx(celnod(ci,3))+nodx(celnod(ci,1)))/2;
                    cx2=(nodx(celnod(ci2,3))+nodx(celnod(ci2,1)))/2;
                    ddxmvy=(xm(m)-cx)/(cx2-cx); 
                    % => Lower right cell > Left cell
                    vy1=S(celvar(ci,2));
                    vy2=S(celvar(ci,3));
                    vy3=(S(celvar(ci2,2))+S(celvar(ci2,3)))/2;
                    vy4=S(celvar(ci2,3));
                end
            end
        else
            % Interpolate Left
            ddxmvy=ddxm+0.5;
            % Left boundary
            if(nodx(celnod(ci,1))==0)
                % No slip
                if(bleft==0)
                    vy1=-S(celvar(ci,2));
                    vy2=-S(celvar(ci,3));
                    vy3=S(celvar(ci,2));
                    vy4=S(celvar(ci,3));
                 else
                    % Free slip
                    vy1=S(celvar(ci,2));
                    vy2=S(celvar(ci,3));
                    vy3=S(celvar(ci,2));
                    vy4=S(celvar(ci,3));
                end
            else
                % Index of the upper left cell
                ci1=nodcel(celnod(ci,1),2);
                % Index of the lower left cell
                ci2=nodcel(celnod(ci,2),1);
                % Upper left cell exist
                if(ci1>0) 
                    % Same size for left and right cell
                    if(celres(ci)==celres(ci1))
                        vy1=S(celvar(ci1,2));
                        vy2=S(celvar(ci1,3));
                        vy3=S(celvar(ci,2));
                        vy4=S(celvar(ci,3));
                    else
                        % Recompute normalized vertical distance
                        cx=(nodx(celnod(ci,3))+nodx(celnod(ci,1)))/2;
                        cx1=(nodx(celnod(ci1,3))+nodx(celnod(ci1,1)))/2;
                        ddxmvy=(xm(m)-cx1)/(cx-cx1); 
                        % Left cell < right cell
                        if(celres(ci)<celres(ci1))
                            % Upper subcell
                            if(ddym<0.5)
                                vy1=S(celvar(ci1,2));
                                vy2=S(celvar(ci1,3));
                                vy3=S(celvar(ci,2));
                                vy4=(S(celvar(ci,2))+S(celvar(ci,3)))/2;
                                ddymvy=ddym*2;
                            else
                                % Lower subcell
                                vy1=S(celvar(ci2,2));
                                vy2=S(celvar(ci2,3));
                                vy3=(S(celvar(ci,2))+S(celvar(ci,3)))/2;
                                vy4=S(celvar(ci,3));
                                ddymvy=(ddym-0.5)*2;
                            end
                        else
                            % Left cell > Right cell
                            % Upper subcell
                            vy1=S(celvar(ci1,2));
                            vy2=(S(celvar(ci1,2))+S(celvar(ci1,3)))/2;
                            vy3=S(celvar(ci,2));
                            vy4=S(celvar(ci,3));
                        end
                    end
                else
                    % Upper left cell does not exist =>
                    % Recompute normalized vertical distance
                    cx=(nodx(celnod(ci,3))+nodx(celnod(ci,1)))/2;
                    cx2=(nodx(celnod(ci2,3))+nodx(celnod(ci2,1)))/2;
                    ddxmvy=(xm(m)-cx2)/(cx-cx2); 
                    % => Lower left cell > Right cell
                    vy1=(S(celvar(ci2,2))+S(celvar(ci2,3)))/2;
                    vy2=S(celvar(ci2,3));
                    vy3=S(celvar(ci,2));
                    vy4=S(celvar(ci,3));
                end
            end
        end
        % Compute vx velocity for the marker
        vym=(1-ddxmvy)*(1-ddymvy)*vy1+...
            (1-ddxmvy)*(ddymvy)*vy2+...
            (ddxmvy)*(1-ddymvy)*vy3+...
            (ddxmvy)*(ddymvy)*vy4;
        
         % Move marker by velocity field
        xm(m)=xm(m)+vxm*timestp;
        ym(m)=ym(m)+vym*timestp;
              
    else
        % Inject new markers
        if(ym(m)<g.ymin)
            ym(m)=ym(m)+bupper(1)*timestp;
        end
    end
end

end
