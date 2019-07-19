% Function varindex_planet()
% This function index variables for agiven cell
% or goes for another level of resolution
% Function return indexed arrays
function[celvar,varnum,varnump]=varindex_planet(ci,varnum,varnump,nodcel,celnod,celvar,celdot)

% Define indexes of unknowns
%       vy2
% vx1  P5,FI   vx4
%       vy3
if(celdot(ci,1)==0)
    % vx1
    % Define cell to the left
    ciL=nodcel(celnod(ci,1),2);
    if(ciL==0)
        ciL=nodcel(celnod(ci,2),1);
    end
    % Add variable    
    if(ciL==0 || celvar(ciL,4)==0)
        % Add new unknown
        varnum=varnum+1;
        celvar(ci,1)=varnum;
        if(ciL>0)
            celvar(ciL,4)=varnum;
        end
    else
        % Use existing variable
        celvar(ci,1)=celvar(ciL,4);
    end
    % vy2
    % Define cell to above
    ciA=nodcel(celnod(ci,1),3);
    if(ciA==0)
        ciA=nodcel(celnod(ci,3),1);
    end
    if(ciA==0 || celvar(ciA,3)==0)
        % Add new unknown
        varnum=varnum+1;
        celvar(ci,2)=varnum;
        if(ciA>0)
            celvar(ciA,3)=varnum;
        end
    else
        % Use existing variable
        celvar(ci,2)=celvar(ciA,3);
    end
    % vy3
    % Define cell to below
    ciB=nodcel(celnod(ci,2),4);
    if(ciB==0)
        ciB=nodcel(celnod(ci,4),2);
    end
    if(ciB==0 || celvar(ciB,2)==0)
        % Add new unknown
        varnum=varnum+1;
        celvar(ci,3)=varnum;
        if(ciB>0)
            celvar(ciB,2)=varnum;
        end
    else
        % Use existing variable
        celvar(ci,3)=celvar(ciB,2);
    end
    % vx4
    % Define cell to the right
    ciR=nodcel(celnod(ci,3),4);
    if(ciR==0)
        ciR=nodcel(celnod(ci,4),3);
    end
    % Add variable    
    if(ciR==0 || celvar(ciR,1)==0)
        % Add new unknown
        varnum=varnum+1;
        celvar(ci,4)=varnum;
        if(ciR>0)
            celvar(ciR,1)=varnum;
        end
    else
        celvar(ci,4)=celvar(ciR,1);
    end
    % P5  
    % Add new unknown
    varnum=varnum+1;
    celvar(ci,5)=varnum;
    % Gravity potential FI 
    % Add new unknown
    varnump=varnump+1;
    celvar(ci,6)=varnump;
else
    % Go to another level of resolution
    for di=1:1:4
        ci1=celdot(ci,di);
        [celvar,varnum,varnump]=varindex_planet(ci1,varnum,varnump,nodcel,celnod,celvar,celdot);
    end
end
end

