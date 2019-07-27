% Function varindex_test2()
% This function index variables for agiven cell
% or goes for another level of resolution
% Function return indexed arrays
function[celvar,varnum]=varindex_test2(ci,varnum,nodcel,celnod,celvar,celdot,celres)

% Define indexes of unknowns
%       vy2
% vx1   P5    vx4
%       vy3
if(celdot(ci,1)==0)
    % vx1
    % Define cell to the left
    ciL=nodcel(celnod(ci,1),2);
    % Add variable    
    if(ciL==0 || celvar(ciL,4)==0 || celres(ciL)~=celres(ci))
        % Add new unknown
        varnum=varnum+1;
        celvar(ci,1)=varnum;
    else
        % Use existing variable
        celvar(ci,1)=celvar(ciL,4);
    end
    % vy2
    % Define cell to above
    ciA=nodcel(celnod(ci,1),3);
    if(ciA==0 || celvar(ciA,3)==0 || celres(ciA)~=celres(ci))
        % Add new unknown
        varnum=varnum+1;
        celvar(ci,2)=varnum;
    else
        % Use existing variable
        celvar(ci,2)=celvar(ciA,3);
    end
    % vy3
    % Define cell to below
    ciB=nodcel(celnod(ci,2),4);
    if(ciB==0 || celvar(ciB,2)==0 || celres(ciB)~=celres(ci))
        % Add new unknown
        varnum=varnum+1;
        celvar(ci,3)=varnum;
    else
        % Use existing variable
        celvar(ci,3)=celvar(ciB,2);
    end
    % vx4
    % Define cell to the right
    ciR=nodcel(celnod(ci,3),4);
    % Add variable    
    if(ciR==0 || celvar(ciR,1)==0  || celres(ciR)~=celres(ci))
        % Add new unknown
        varnum=varnum+1;
        celvar(ci,4)=varnum;
    else
        celvar(ci,4)=celvar(ciR,1);
    end
    % P5  
    % Add new unknown
    varnum=varnum+1;
    celvar(ci,5)=varnum;
else
    % Go to another level of resolution
    for di=1:1:4
        ci1=celdot(ci,di);
        [celvar,varnum]=varindex_test2(ci1,varnum,nodcel,celnod,celvar,celdot,celres);
    end
end
end

