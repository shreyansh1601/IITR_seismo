% Computation and visualisation of viscosity profile across the lithosphere

% Clearing all variables and arrays
clear all;
% Clearing figures
clf;

% Flow law parameters
AD=2.5e-17; %1/Pa^n/s, 
n=3.5;
Ea=532000; %J/mol 
R=8.314; %J/mol/K

% Define correction coefficient F2 
% for a strain rate based viscosity formulation
F2=1/2^((n-1)/n)/3^((n+1)/2/n)

% Strain Rate
eps=1e-14; %1/s


% Defining temperature profile
L=100000; % Length of the profile, m
ttop=400; % T o the top, oC
tbottom=1200; % T at the bottom, oC
ynum=101; % Number of points in the profile
y=0:L/(ynum-1):L; % depth vector, m
T=ttop:(tbottom-ttop)/(ynum-1):tbottom; % Temperature vector, oC

% Compute viscosity profile 
eta=ones(ynum,1); % viscosity array
for i=1:1:ynum
    % Compute and check activation energy exponent
    eterm=Ea/n/R/(T(i)+273);
    if (eterm>100) 
        eterm=100; 
    end
    eterm=exp(eterm);
    % Compute viscosity
    eta(i)=F2/AD^(1/n)/eps^((n-1)/n)*eterm;
end


% Defining new figure
figure(1)

% Plotting vertical  temperature profile
subplot(1,2,1);     % defining 1st plotting area in the 3x1 figure
plot(T,y/1000,'k'); % making plot: 'k' = black color
box on;             % making box around the plot
title('Temperature profile');   % title of the plot
xlabel('T, oC');        % title of horizontal axis
ylabel('Depth, km'); % title of vertical axis
axis ij;
% Plotting vertical log viscosity profile
subplot(1,2,2);     % defining 2nd plotting area in the 3x1 figure
plot(log10(eta),y/1000,'k'); % making plot: 'k' = black color
box on;             % making box around the plot
title('log viscosity profile');   % title of the plot
xlabel('ETA, Pa s');        % title of horizontal axis
ylabel('Depth, km'); % title of vertical axis
axis ij;

