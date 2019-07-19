% Computation of visco-elasto-plastic stress buildup 
% and associated viscous, elastic and plastic strain rate evolution with time

% Clear all arrays
clear all; 

% Clear all figures
clf;

tnum=100; %number of timesteps
tstp=1e+3*(365.25*24*3600); % timestep, s

% Time vector
t=0:tstp:(tstp*(tnum-1));

% Stress vector 
sij=zeros(1,tnum);
% Strain rate vectors 
eijviscous=zeros(1,tnum); % viscous strain rate
eijelastic=zeros(1,tnum); % elastic strain rate
eijplastic=zeros(1,tnum); % plastic strain rate
% Parameters for 6 cases
% Initial stress, Pa
s0ij  = 0;
% Strain rate, 1/s
eij   = 1e-14;
% Viscosity, Pa s
eta   = 1e+22;
% Shear modulus, Pa
mu    = 1e+10;
% Yeld stress, Pa
syield = 1.5e+8;

% Compute stress(time) curves
% Time cycle 
for k=1:1:tnum
    % Stress, Pa
    sij(k)=2*eta*eij+(s0ij-2*eta*eij)*exp(-t(k)*mu/eta);
    % Viscoplastic deformation (on the yield stress)
    if (sij(k)>syield)
        sij(k)=syield;
        eijviscous(k)=sij(k)/2/eta; % viscous strain rate
        eijelastic(k)=0; % elastic strain rate =0 since stress does not change with time (Sij=Syield=const)
        eijplastic(k)=eij-eijviscous(k); % plastic strain rate (since elastic strain rate =0 on the yield)
    else
    % Viscoelastic deformation (below the yield stress)
        eijviscous(k)=sij(k)/2/eta; % viscous strain rate
        eijelastic(k)=eij-eijviscous(k); % elastic strain rate (since plastic strain rate = 0 below yield stress)
        eijplastic(k)=0; % plastic strain rate = 0 below yeild
    end
end

% Recompute time to years
t=t/(365.25*24*3600);
% Recompute stress to MPa
sij=sij*1e-6;
% Plot results 
figure(1);
% Plot stress
subplot(1,2,1);
plot(t,sij,'r');
% Plot strain rates
subplot(1,2,2);
% Viscous
plot(t,eijviscous,'k');
hold on;
% Elastic
plot(t,eijelastic,'r');
% Plastic
plot(t,eijplastic,'b');
hold off;