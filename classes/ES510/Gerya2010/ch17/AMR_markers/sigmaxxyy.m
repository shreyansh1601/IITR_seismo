% Function sigmaxxyy_test()
% This function compose EPSILONxx, EPSILONyy, SIGMA'xx, SIGMA'yy for specific cell
% Function return modified matrix
function[exx1,sxx1,eyy1,syy1]=sigmaxxyy(ci,nodx,nody,celnod,celvar,celeta,S)

%       vy2
%       P5 
%       vy3
% Cell x-size
dx=nodx(celnod(ci,3))-nodx(celnod(ci,1));
dy=nody(celnod(ci,2))-nody(celnod(ci,1));
% Compute EPSILONxx, EPSILONyy, SIGMA'xx, SIGMA'yy 
exx1=(S(celvar(ci,4))-S(celvar(ci,1)))/dx; % EPSILONxx
eyy1=(S(celvar(ci,3))-S(celvar(ci,2)))/dy; % EPSILONyy
sxx1=2*celeta(ci)*exx1; % SIGMA'xx
syy1=2*celeta(ci)*eyy1; % SIGMA'yy
end

