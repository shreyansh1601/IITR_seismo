% Computation and visualisation of viscosity map 
% in temperature - log stress coordinates

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
% for a stress based viscosity formulation
F1=1/3^((n+1)/2)


% Defining temperature variations
tmin=400;  % minimal T, oC
tmax=1400; % maximal T, oC
tnum=51; % Number of points along T axis
T=tmin:(tmax-tmin)/(tnum-1):tmax; % Temperature vector, oC
% Defining temperature variations
smin=3;   % minimal log of stress, Pa
smax=9;   % maximal log of stress, Pa
snum=51; % Number of points along log stress axis
S=smin:(smax-smin)/(snum-1):smax; % log stress vector, Pa

% Compute viscosity profile 
eta=ones(snum,tnum); % viscosity array
for i=1:1:snum
    for j=1:1:tnum
        % Compute and check activation energy exponent
        eterm=Ea/R/(T(j)+273);
        if (eterm>100) 
            eterm=100; 
        end
        eterm=exp(eterm);
        % Compute viscosity
        eta(i,j)=F1/AD/(10^S(i))^(n-1)*eterm;
    end
end


% Defining new figure
figure(1)

% Plotting viscosity map
pcolor(T,S,log10(eta));
shading interp;
colorbar;
box on;             
title('log viscosity map'); 
xlabel('T, ^oC');        
ylabel('log stress, Pa'); 


