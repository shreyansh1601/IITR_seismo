% Function celinter()
% This function interpolate a property 
% from cell centers for given location 
% Function return the interpolated value
function[ival]=celinter(ci,cval,xx,yy,nodx,nody,nodcel,celnod,btleft,btright,btupper,btlower,xsize,ysize)

% Arrays
cii=zeros(13,1);
civ=zeros(13,1);
cix=zeros(13,1);
ciy=zeros(13,1);
cidx=zeros(13,1);
cidy=zeros(13,1);
%  10|  6  |  7  | 12
% ---+-----+-----+----
%  2 |  o        | 4
% ---+     1 ci  +---
%  3 |           | 5
% ---+-----+-----+----
%  11|  8  |  9  | 13
%
% Process all cells

% Cell 1 
% Index
cii(1)=ci;
% Cell coordinates
cix(1)=(nodx(celnod(ci,3))+nodx(celnod(ci,1)))/2;
ciy(1)=(nody(celnod(ci,2))+nody(celnod(ci,1)))/2;
% Cell size
cidx(1)=nodx(celnod(ci,3))-nodx(celnod(ci,1));
cidy(1)=nody(celnod(ci,2))-nody(celnod(ci,1));
% Value
civ(1)=cval(ci);
 
% Cell 2
if(nodx(celnod(ci,1))==0)
    % Index
    cii(2)=-1;
    % Cell coordinates
    cix(2)=cix(1)-cidx(1);
    ciy(2)=ciy(1);
    % Cell size
    cidx(2)=cidx(1);
    cidy(2)=cidy(1);
    % Value
    % Boundary condition
    if(btleft<0)
        % Horizontal symmetry
        civ(2)=civ(1);
    else
        % Constant value
        civ(2)=2*btleft-civ(1);
    end
else
   % Cell index
   ci2=nodcel(celnod(ci,1),2);
   if(ci2>0)
        % Index
        cii(2)=ci2;
        % Cell coordinates
        cix(2)=(nodx(celnod(ci2,3))+nodx(celnod(ci2,1)))/2;
        ciy(2)=(nody(celnod(ci2,2))+nody(celnod(ci2,1)))/2;
        % Cell size
        cidx(2)=nodx(celnod(ci2,3))-nodx(celnod(ci2,1));
        cidy(2)=nody(celnod(ci2,2))-nody(celnod(ci2,1));
        % Value
        civ(2)=cval(ci2);
   end
   %
   % Cell 3
   % Cell index
   ci3=nodcel(celnod(ci,2),1);
   if(ci3>0 && ci3~=cii(2))
        % Index
        cii(3)=ci3;
        % Cell coordinates
        cix(3)=(nodx(celnod(ci3,3))+nodx(celnod(ci3,1)))/2;
        ciy(3)=(nody(celnod(ci3,2))+nody(celnod(ci3,1)))/2;
        % Cell size
        cidx(3)=nodx(celnod(ci3,3))-nodx(celnod(ci3,1));
        cidy(3)=nody(celnod(ci3,2))-nody(celnod(ci3,1));
        % Value
        civ(3)=cval(ci3);
   end
end



% Cell 4
if(nodx(celnod(ci,3))==xsize)
    % Index
    cii(4)=-1;
    % Cell coordinates
    cix(4)=cix(1)+cidx(1);
    ciy(4)=ciy(1);
    % Cell size
    cidx(4)=cidx(1);
    cidy(4)=cidy(1);
    % Value
    % Boundary condition
    if(btright<0)
        % Horizontal symmetry
        civ(4)=civ(1);
    else
        % Constant value
        civ(4)=2*btright-civ(1);
    end
else
   % Cell index
   ci4=nodcel(celnod(ci,3),4);
   if(ci4>0)
        % Index
        cii(4)=ci4;
        % Cell coordinates
        cix(4)=(nodx(celnod(ci4,3))+nodx(celnod(ci4,1)))/2;
        ciy(4)=(nody(celnod(ci4,2))+nody(celnod(ci4,1)))/2;
        % Cell size
        cidx(4)=nodx(celnod(ci4,3))-nodx(celnod(ci4,1));
        cidy(4)=nody(celnod(ci4,2))-nody(celnod(ci4,1));
        % Value
        civ(4)=cval(ci4);
   end
   %
   % Cell 5
   % Cell index
   ci5=nodcel(celnod(ci,4),3);
   if(ci5>0 && ci5~=cii(4))
        % Index
        cii(5)=ci5;
        % Cell coordinates
        cix(5)=(nodx(celnod(ci5,3))+nodx(celnod(ci5,1)))/2;
        ciy(5)=(nody(celnod(ci5,2))+nody(celnod(ci5,1)))/2;
        % Cell size
        cidx(5)=nodx(celnod(ci5,3))-nodx(celnod(ci5,1));
        cidy(5)=nody(celnod(ci5,2))-nody(celnod(ci5,1));
        % Value
        civ(5)=cval(ci5);
   end
end


% Cell 6
if(nody(celnod(ci,1))==0)
    % Index
    cii(6)=-1;
    % Cell coordinates
    cix(6)=cix(1);
    ciy(6)=ciy(1)-cidy(1);
    % Cell size
    cidx(6)=cidx(1);
    cidy(6)=cidy(1);
    % Value
    % Boundary condition
    if(btupper<0)
        % Horizontal symmetry
        civ(6)=civ(1);
    else
        % Constant value
        civ(6)=2*btupper-civ(1);
    end
else
   % Cell index
   ci6=nodcel(celnod(ci,1),3);
   if(ci6>0)
        % Index
        cii(6)=ci6;
        % Cell coordinates
        cix(6)=(nodx(celnod(ci6,3))+nodx(celnod(ci6,1)))/2;
        ciy(6)=(nody(celnod(ci6,2))+nody(celnod(ci6,1)))/2;
        % Cell size
        cidx(6)=nodx(celnod(ci6,3))-nodx(celnod(ci6,1));
        cidy(6)=nody(celnod(ci6,2))-nody(celnod(ci6,1));
        % Value
        civ(6)=cval(ci6);
   end
   %
   % Cell 7
   % Cell index
   ci7=nodcel(celnod(ci,3),1);
   if(ci7>0 && ci7~=cii(6))
        % Index
        cii(7)=ci7;
        % Cell coordinates
        cix(7)=(nodx(celnod(ci7,3))+nodx(celnod(ci7,1)))/2;
        ciy(7)=(nody(celnod(ci7,2))+nody(celnod(ci7,1)))/2;
        % Cell size
        cidx(7)=nodx(celnod(ci7,3))-nodx(celnod(ci7,1));
        cidy(7)=nody(celnod(ci7,2))-nody(celnod(ci7,1));
        % Value
        civ(7)=cval(ci7);
   end
end

% Cell 8
if(nody(celnod(ci,2))==ysize)
    % Index
    cii(8)=-1;
    % Cell coordinates
    cix(8)=cix(1);
    ciy(8)=ciy(1)+cidy(1);
    % Cell size
    cidx(8)=cidx(1);
    cidy(8)=cidy(1);
    % Value
    % Boundary condition
    if(btlower<0)
        % Horizontal symmetry
        civ(8)=civ(1);
    else
        % Constant value
        civ(8)=2*btlower-civ(1);
    end
else
   % Cell index
   ci8=nodcel(celnod(ci,2),4);
   if(ci8>0)
        % Index
        cii(8)=ci8;
        % Cell coordinates
        cix(8)=(nodx(celnod(ci8,3))+nodx(celnod(ci8,1)))/2;
        ciy(8)=(nody(celnod(ci8,2))+nody(celnod(ci8,1)))/2;
        % Cell size
        cidx(8)=nodx(celnod(ci8,3))-nodx(celnod(ci8,1));
        cidy(8)=nody(celnod(ci8,2))-nody(celnod(ci8,1));
        % Value
        civ(8)=cval(ci8);
   end
   %
   % Cell 9
   % Cell index
   ci9=nodcel(celnod(ci,4),2);
   if(ci9>0 && ci9~=cii(8))
        % Index
        cii(9)=ci9;
        % Cell coordinates
        cix(9)=(nodx(celnod(ci9,3))+nodx(celnod(ci9,1)))/2;
        ciy(9)=(nody(celnod(ci9,2))+nody(celnod(ci9,1)))/2;
        % Cell size
        cidx(9)=nodx(celnod(ci9,3))-nodx(celnod(ci9,1));
        cidy(9)=nody(celnod(ci9,2))-nody(celnod(ci9,1));
        % Value
        civ(9)=cval(ci9);
   end
end

% Cell 10
if(nodx(celnod(ci,1))==0)
    % Index
    cii(10)=-1;
    % Cell coordinates
    cix(10)=cix(6)-cidx(6);
    ciy(10)=ciy(6);
    % Cell size
    cidx(10)=cidx(6);
    cidy(10)=cidy(6);
    % Value
    % Boundary condition
    if(btleft<0)
        % Horizontal symmetry
        civ(10)=civ(6);
    else
        % Constant value
        civ(10)=2*btleft-civ(6);
    end
else
   % Cell index
   ci10=nodcel(celnod(ci,1),1);
   if(ci10>0)
        % Index
        cii(10)=ci10;
        % Cell coordinates
        cix(10)=(nodx(celnod(ci10,3))+nodx(celnod(ci10,1)))/2;
        ciy(10)=(nody(celnod(ci10,2))+nody(celnod(ci10,1)))/2;
        % Cell size
        cidx(10)=nodx(celnod(ci10,3))-nodx(celnod(ci10,1));
        cidy(10)=nody(celnod(ci10,2))-nody(celnod(ci10,1));
        % Value
        civ(10)=cval(ci10);
   end
end


% Cell 11
if(nodx(celnod(ci,1))==0)
    % Index
    cii(11)=-1;
    % Cell coordinates
    cix(11)=cix(8)-cidx(8);
    ciy(11)=ciy(8);
    % Cell size
    cidx(11)=cidx(8);
    cidy(11)=cidy(8);
    % Value
    % Boundary condition
    if(btleft<0)
        % Horizontal symmetry
        civ(11)=civ(8);
    else
        % Constant value
        civ(11)=2*btleft-civ(8);
    end
else
   % Cell index
   ci11=nodcel(celnod(ci,2),2);
   if(ci11>0)
        % Index
        cii(11)=ci11;
        % Cell coordinates
        cix(11)=(nodx(celnod(ci11,3))+nodx(celnod(ci11,1)))/2;
        ciy(11)=(nody(celnod(ci11,2))+nody(celnod(ci11,1)))/2;
        % Cell size
        cidx(11)=nodx(celnod(ci11,3))-nodx(celnod(ci11,1));
        cidy(11)=nody(celnod(ci11,2))-nody(celnod(ci11,1));
        % Value
        civ(11)=cval(ci11);
   end
end


% Cell 12
if(nodx(celnod(ci,3))==xsize)
    ciR=7;
    if(cii(7)==0)
        ciR=6;
    end
    % Index
    cii(12)=-1;
    % Cell coordinates
    cix(12)=cix(ciR)+cidx(ciR);
    ciy(12)=ciy(ciR);
    % Cell size
    cidx(12)=cidx(ciR);
    cidy(12)=cidy(ciR);
    % Value
    % Boundary condition
    if(btright<0)
        % Horizontal symmetry
        civ(12)=civ(ciR);
    else
        % Constant value
        civ(12)=2*btright-civ(ciR);
    end
else
   % Cell index
   ci12=nodcel(celnod(ci,3),3);
   if(ci12>0)
        % Index
        cii(12)=ci12;
        % Cell coordinates
        cix(12)=(nodx(celnod(ci12,3))+nodx(celnod(ci12,1)))/2;
        ciy(12)=(nody(celnod(ci12,2))+nody(celnod(ci12,1)))/2;
        % Cell size
        cidx(12)=nodx(celnod(ci12,3))-nodx(celnod(ci12,1));
        cidy(12)=nody(celnod(ci12,2))-nody(celnod(ci12,1));
        % Value
        civ(12)=cval(ci12);
   end
end

% Cell 13
if(nodx(celnod(ci,3))==xsize)
    ciR=9;
    if(cii(9)<0)
        ciR=8;
    end
    % Index
    cii(13)=-1;
    % Cell coordinates
    cix(13)=cix(ciR)+cidx(ciR);
    ciy(13)=ciy(ciR);
    % Cell size
    cidx(13)=cidx(ciR);
    cidy(13)=cidy(ciR);
    % Value
    % Boundary condition
    if(btright==0)
        % Horizontal symmetry
        civ(13)=civ(ciR);
    else
        % Constant value
        civ(13)=2*btright-civ(ciR);
    end
else
   % Cell index
   ci13=nodcel(celnod(ci,4),4);
   if(ci13>0)
        % Index
        cii(13)=ci13;
        % Cell coordinates
        cix(13)=(nodx(celnod(ci13,3))+nodx(celnod(ci13,1)))/2;
        ciy(13)=(nody(celnod(ci13,2))+nody(celnod(ci13,1)))/2;
        % Cell size
        cidx(13)=nodx(celnod(ci13,3))-nodx(celnod(ci13,1));
        cidy(13)=nody(celnod(ci13,2))-nody(celnod(ci13,1));
        % Value
        civ(13)=cval(ci13);
   end
end

% Compute interpolated value
ival=0;
wt=0;
% Update the value from 13 cells
for i=1:1:13
    if(cii(i)~=0)
        % Normalized distances
        dxx=1-abs(xx-cix(i))/cidx(i);
        dyy=1-abs(yy-ciy(i))/cidy(i);
        % Check distances, Update the value
        if(dxx>0 && dyy>0)
            ival=ival+civ(i)*dxx*dyy;
            wt=wt+dxx*dyy;
        end
    end
end
% Compute interpolated value
ival=ival/wt;
end
