% Solution of 1D Poisson equation d2FI/dx=1 
% with finite differences on a regular grid
% using direct solver ‘\’

% Clean all variables
clear all;
% Clear all figures
clf;

% Numerical model parameters
xsize=1000000; % Model size, m
xnum=1000;    % Number of nodes
xstp=xsize/(xnum-1); % Grid step

% Right part of Poisson equation is = 1
% Boundary conditions: FI=0

% Matrix of coefficients initialization
L=sparse(xnum,xnum);
% Vector of right part initialization
R=zeros(xnum,1);

% Composing matrix of coefficients L()
% and vector (column) of right parts R()
% First point: FI=0
L(1,1)=1;
R(1,1)=0;
% Intermediate points
for i=2:1:xnum-1
    % d2FI/dx2=1
    % (FI(i-1)-2*FI(i)+FI(i+1))/dx^2=1
    L(i,i-1)=1/xstp^2;
    L(i,i)=-2/xstp^2;
    L(i,i+1)=1/xstp^2;
    R(i,1)=1;
end
% Last point: FI=0
L(xnum,xnum)=1;
R(xnum,1)=0;

% Obtaining vector (line) of solutions S()
S=L\R;

%Creating vector for nodal point positions
x=0:xstp:xsize;

% Open figure
figure(7);
% Plotting solutions
plot(x/1000,S);











