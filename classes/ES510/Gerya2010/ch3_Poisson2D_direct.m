% Solution of 2D Poisson equation d2FI/dx2+d2FI/dy2=1 
% with finite differences on a regular grid
% using direct solver ‘\’

% Clean all variables
clear all;
% Clear all figures
clf;

% Numerical model parameters
% Model size, m
xsize=1000; % Horizontal
ysize=1500; % Vertical

% Numbers of nodes
xnum=31; % Horizontal
ynum=41; % Vertical
% Grid step
xstp=xsize/(xnum-1); % Horizontal
ystp=ysize/(ynum-1); % Vertical

% Matrix of coefficients initialization
L=sparse(xnum*ynum,xnum*ynum);
% Vector of right part initialization
R=zeros(xnum*ynum,1);

%Composing matrix of coefficients L()
%and vector (column) of right parts R()
% Boundary conditions: FI=0
% Process all Grid points
for i=1:1:ynum
  for j=1:1:xnum
    % Global index for current node
    k=(j-1)*ynum+i;
    %Boundary nodes
    if(i==1 || i==ynum || j==1 || j==xnum)
        L(k,k)=1;
        R(k,1)=0;
    else
    %Internal nodes
        %d2FI/dx2+d2FI/dy2=1
        L(k,k-ynum)=1/xstp^2;
        L(k,k-1)=1/ystp^2;
        L(k,k)=-2/xstp^2-2/ystp^2;
        L(k,k+1)=1/ystp^2;
        L(k,k+ynum)=1/xstp^2;
        R(k,1)=1;
    end
  end
end
    
%Obtaining vector of solutions S()
S=L\R;

% Reload solutions S() to 2D array FI()
FI=zeros(ynum,xnum);
% Process all Grid points
for i=1:1:ynum
  for j=1:1:xnum
    % Global index for current node
    k=(j-1)*ynum+i;
    FI(i,j)=S(k);
  end
end

% Making vectors for nodal points positions
x=0:xstp:xsize; % Horizontal
y=0:ystp:ysize; % Vertical

% Making new figure
figure(7);
%Plotting solution
%pcolor(PHI);
surf(x,y,FI);
light;
shading interp;
colorbar;
lighting phong;
