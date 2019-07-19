% Function sigmaxx()
% This function compose SIGMAxx=SIGMA'xx-P for specific cell
% and add it to the global matrix
% Function return modified matrix
function[L]=sigmaxx(kk,ci,nodx,celnod,celvar,celeta,vkf,pkf,L)

% vx1   P5   vx4
% Cell x-size
dx=nodx(celnod(ci,3))-nodx(celnod(ci,1));
% Matrix update
L(kk,celvar(ci,1))=L(kk,celvar(ci,1))-vkf*2*celeta(ci)/dx; % vx1
L(kk,celvar(ci,4))=L(kk,celvar(ci,4))+vkf*2*celeta(ci)/dx; % vx4
L(kk,celvar(ci,5))=L(kk,celvar(ci,5))-pkf; % P5         
end

