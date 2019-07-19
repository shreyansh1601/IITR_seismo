% Function mergecelbest()
% This function split given parent cell ci
% to 4 subcells and add additional nodes
% Function return modified nodal and cell arrays
function [celcor,celdot,celnod,nodcel,recnum,celrec,rennum,nodrec]=mergecelbest(ci,celdot,celnod,nodcel,celcor,recnum,celrec,rennum,nodrec)


% Mark cells for recycling
recnum=recnum+1;
celrec(recnum)=celdot(ci,1);
% 4 nodes connected to the cell
% 1    3
%   ci
% 2    4
% Delete daughter cells from nodes
for cid=1:1:4
    for ni=1:1:4
        if(ni~=cid)
            % Remove daughter cells from nodes
            nodcel(celnod(celdot(ci,cid),ni),5-ni)=0;
            % Recycle empty nodes
            recyn=1;
            for nic=1:1:4
                if(nodcel(celnod(celdot(ci,cid),ni),nic)>0)
                    recyn=0;
                end
            end
            % Node is not connected to any cell
            if(recyn==1)
                rennum=rennum+1;
                nodrec(rennum)=celnod(celdot(ci,cid),ni);
            end
        else
            % Add Parent Cell to nodes
            nodcel(celnod(ci,cid),5-cid)=ci;
        end
    end
    % Remove coarsening mark
    celcor(celdot(ci,cid))=0;
    % Set non-working cell mark
    celdot(celdot(ci,cid),1)=-1;
end
% Delete subcel reference
celdot(ci,1)=0;
end

