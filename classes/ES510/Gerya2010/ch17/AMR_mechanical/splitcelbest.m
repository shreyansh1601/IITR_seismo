% Function splitcelbest()
% This function split given cell ci
% to 4 subcells and add additional nodes
% using recycled cells if available
% Function return modified nodal and cell arrays
function[celnum,nodnum,celdot,celres,celpar,celnod,nodcel,nodx,nody,recnum,celrec,rennum,nodrec]=...
     splitcelbest(ci,celnum,nodnum,celdot,celres,celpar,celnod,nodcel,nodx,nody,recnum,celrec,rennum,nodrec)

 % Resolution level for daughter cells
reslev1=celres(ci)+1;

% Cell center coordinates
cx=(nodx(celnod(ci,1))+nodx(celnod(ci,4)))/2;
cy=(nody(celnod(ci,1))+nody(celnod(ci,4)))/2;

% Save cell number
celnum0=celnum;
% Use old daughter cells
newcel=1;
if(recnum>0)
    newcel=0;
    celnum=celrec(recnum)-1;
    recnum=recnum-1;
end

% Add 4 doughter cells
% 1  3
% 2  4
for cid=1:1:4
    % Daughter cells
    celdot(ci,cid)=celnum+cid;
    % Resolution level
    celres(celnum+cid)=reslev1;
    % Working cell mark
    celdot(celnum+cid,1)=0;
    % Parent cell
    celpar(celnum+cid,1)=ci;
end

% Create central node in the cell
if(rennum==0)
    nodnum=nodnum+1;
    nodcur=nodnum;
else
    nodcur=nodrec(rennum);
    rennum=rennum-1;
end
nodx(nodcur)=cx; % x coordinate
nody(nodcur)=cy; % y coordinate
% Add  parent and central nodes
for cid=1:1:4
    % Parent node
    celnod(celnum+cid,cid)=celnod(ci,cid);
    nodcel(celnod(celnum+cid,cid),5-cid)=celnum+cid;
    % Central node
    celnod(celnum+cid,5-cid)=nodcur;
    nodcel(nodcur,cid)=celnum+cid;
end

% Add remaining nodes
% Cell 1--------------------------------------
% 4 nodes connected to the cell
%      3
%   ci
% 2    
% Node 2
if(nodcel(celnod(ci,1),2)>0 && celres(nodcel(celnod(ci,1),2))==reslev1)
    celnod(celnum+1,2)=celnod(nodcel(celnod(ci,1),2),4);
else
    if(rennum==0)
        nodnum=nodnum+1;
        nodcur=nodnum;
    else
        nodcur=nodrec(rennum);
        rennum=rennum-1;
    end
    nodx(nodcur)=nodx(celnod(ci,1)); % x coordinate
    nody(nodcur)=cy; % y coordinate
    celnod(celnum+1,2)=nodcur;
end
nodcel(celnod(celnum+1,2),3)=celnum+1;
% Node 3
if(nodcel(celnod(ci,1),3)>0 && celres(nodcel(celnod(ci,1),3))==reslev1)
    celnod(celnum+1,3)=celnod(nodcel(celnod(ci,1),3),4);
else
    if(rennum==0)
        nodnum=nodnum+1;
        nodcur=nodnum;
    else
        nodcur=nodrec(rennum);
        rennum=rennum-1;
    end
    nodx(nodcur)=cx; % x coordinate
    nody(nodcur)=nody(celnod(ci,1)); % y coordinate
    celnod(celnum+1,3)=nodcur;
end
nodcel(celnod(celnum+1,3),2)=celnum+1;

% Cell 2--------------------------------------
% 4 nodes connected to the cell
% 1    
%   ci
%      4
% Node 1
celnod(celnum+2,1)=celnod(celnum+1,2);
nodcel(celnod(celnum+2,1),4)=celnum+2;
% Node 4
if(nodcel(celnod(ci,2),4)>0 && celres(nodcel(celnod(ci,2),4))==reslev1)
    celnod(celnum+2,4)=celnod(nodcel(celnod(ci,2),4),3);
else
    if(rennum==0)
        nodnum=nodnum+1;
        nodcur=nodnum;
    else
        nodcur=nodrec(rennum);
        rennum=rennum-1;
    end
    nodx(nodcur)=cx; % x coordinate
    nody(nodcur)=nody(celnod(ci,2)); % y coordinate
    celnod(celnum+2,4)=nodcur;
end
nodcel(celnod(celnum+2,4),1)=celnum+2;

% Cell 3--------------------------------------
% 4 nodes connected to the cell
% 1     
%   ci
%      4
% Node 1
celnod(celnum+3,1)=celnod(celnum+1,3);
nodcel(celnod(celnum+3,1),4)=celnum+3;
% Node 4
if(nodcel(celnod(ci,3),4)>0 && celres(nodcel(celnod(ci,3),4))==reslev1)
    celnod(celnum+3,4)=celnod(nodcel(celnod(ci,3),4),2);
else
    if(rennum==0)
        nodnum=nodnum+1;
        nodcur=nodnum;
    else
        nodcur=nodrec(rennum);
        rennum=rennum-1;
    end
    nodx(nodcur)=nodx(celnod(ci,3)); % x coordinate
    nody(nodcur)=cy; % y coordinate
    celnod(celnum+3,4)=nodcur;
end
nodcel(celnod(celnum+3,4),1)=celnum+3;

% Cell 4--------------------------------------
% 4 nodes connected to the cell
%      3
%   ci
% 2    
% Node 2
celnod(celnum+4,2)=celnod(celnum+2,4);
nodcel(celnod(celnum+4,2),3)=celnum+4;
% Node 3
celnod(celnum+4,3)=celnod(celnum+3,4);
nodcel(celnod(celnum+4,3),2)=celnum+4;
% ----------------------------
% Add number of cells
if(newcel==1)
    celnum=celnum+4;
else
    celnum=celnum0;
end
end

