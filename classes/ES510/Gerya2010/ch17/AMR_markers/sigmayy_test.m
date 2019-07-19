% Function sigmayy_test()
% This function compose SIGMAyy=SIGMA'yy-P for specific cell
% and add it to the global matrix
% Function return modified matrix
function[L]=sigmayy_test(kk,ci,nody,celnod,celvar,celeta,vkf,pkf,L)

%       vy2
%       P5 
%       vy3
% Cell x-size
dy=nody(celnod(ci,2))-nody(celnod(ci,1));
% Matrix update
L(kk,celvar(ci,2))=L(kk,celvar(ci,2))-vkf*2*celeta(ci)/dy; % vy2
L(kk,celvar(ci,3))=L(kk,celvar(ci,3))+vkf*2*celeta(ci)/dy; % vy3
L(kk,celvar(ci,5))=L(kk,celvar(ci,5))-pkf; % P5         
end

