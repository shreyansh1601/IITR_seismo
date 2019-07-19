% Computation and visualisation of viscosity maps in temperature - log stress coordinates
% for a combination of dislocation and diffusion creep;
% comparison of wet and dry olivine rheology (Karato and Wu, 1993)

% Clearing all variables and arrays
clear all;
% Clearing figures
clf;

% Shear Modulus, Pa
mu=8e+10; %80 GPa
% Burgers vector, m
b=5e-10; % 0.5 nm
% Grain size, m
h=1e-3; % 1 mm
% Dislocation creep	[dry wet]	
A_dis= [3.5e+22	2.0e+18]; % 1/s	
n_dis= [3.5	    3      ]; % stress exponent
Ea_dis=[540000	430000 ]; % J/mol
% Diffusion creep		
A_dif= [8.7e+15	5.3e+15]; %1/s
m_dif= [-2.5	-2.5   ]; % grain size exponent
Ea_dif=[300000	240000 ]; % J/mol


R=8.314; %Gas constant, J/mol/K

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
etadry=ones(snum,tnum); % wet viscosity array
etawet=ones(snum,tnum); % wet viscosity array
for i=1:1:snum
    for j=1:1:tnum
        % Dry rheology
        % Dry dislocation creep
        % Compute and check activation energy exponent
        eterm=Ea_dis(1)/R/(T(j)+273);
        if (eterm>100) 
            eterm=100; 
        end
        eterm=exp(-eterm);
        % Compute dry dislocation creep viscosity
        eta_dis=10^S(i)/(2*A_dis(1)*(10^S(i)/mu)^n_dis(1)*eterm);
        % Dry diffusion creep
        % Compute and check activation energy exponent
        eterm=Ea_dif(1)/R/(T(j)+273);
        if (eterm>100) 
            eterm=100; 
        end
        eterm=exp(-eterm);
        % Compute dry diffusion creep viscosity
        eta_dif=10^S(i)/(2*A_dif(1)*(h/b)^m_dif(1)*(10^S(i)/mu)*eterm);
        % Compute effective viscosity for dry creep
        etadry(i,j)=1/(1/eta_dis+1/eta_dif);
        % Wet rheology
        % Wet dislocation creep
        % Compute and check activation energy exponent
        eterm=Ea_dis(2)/R/(T(j)+273);
        if (eterm>100) 
            eterm=100; 
        end
        eterm=exp(-eterm);
        % Compute wet dislocation creep viscosity
        eta_dis=10^S(i)/(2*A_dis(2)*(10^S(i)/mu)^n_dis(2)*eterm);
        % Wet diffusion creep
        % Compute and check activation energy exponent
        eterm=Ea_dif(2)/R/(T(j)+273);
        if (eterm>100) 
            eterm=100; 
        end
        eterm=exp(-eterm);
        % Compute wet diffusion creep viscosity
        eta_dif=10^S(i)/(2*A_dif(2)*(h/b)^m_dif(2)*(10^S(i)/mu)*eterm);
        % Compute effective viscosity for dry creep
        etawet(i,j)=1/(1/eta_dis+1/eta_dif);
    end
end


% Defining new figure
figure(1)

% Plotting dry viscosity map
subplot(1,3,1)
pcolor(T,S,log10(etadry));
shading interp;
colorbar;
box on;             
title('log viscosity map (dry)'); 
xlabel('T, ^oC');        
ylabel('log stress, Pa'); 
% Plotting wet viscosity map
subplot(1,3,2)
pcolor(T,S,log10(etawet));
shading interp;
colorbar;
box on;             
title('log viscosity map (wet)'); 
xlabel('T, ^oC');        
ylabel('log stress, Pa'); 
% Plotting dry/wet viscosity ratio
subplot(1,3,3)
pcolor(T,S,log10(etadry./etawet));
shading interp;
colorbar;
box on;             
title('log viscosity contrast map (dry/wet)'); 
xlabel('T, ^oC');        
ylabel('log stress, Pa'); 




