% Computation of visco-elastic stress buildup/relaxation with time

% Clear all arrays
clear all; 

% Clear all figures
clf;

tnum=100; %number of timesteps
tstp=3e+2*(365.25*24*3600); % timestep, s

% Time vector
t=0:tstp:(tstp*(tnum-1));

% Stress vector for 6 cases
sij=zeros(6,tnum);
% Parameters for 6 cases
%       (1)   (2)   (3)   (4)   (5)   (6)
% Initial stress, Pa
s0ij = [0     1e+8  0     0     0     0    ];
% Strain rate, 1/s
eij  = [1e-14 1e-14 1e-15 1e-14 1e-14 1e-14];
% Viscosity, Pa s
eta  = [1e+21 1e+21 1e+21 1e+22 1e+21 1e+22];
% Shear modulus, Pa
mu   = [1e+10 1e+10 1e+10 1e+10 1e+11 1e+11];
% Compute stress(time) curves
% Time cycle 
for k=1:1:tnum
    % Cases
    for n=1:1:6
        % Stress, Pa
        sij(n,k)=2*eta(n)*eij(n)+(s0ij(n)-2*eta(n)*eij(n))*exp(-t(k)*mu(n)/eta(n));
    end
end

% Recompute time to years
t=t/(365.25*24*3600);
% Recompute stress to MPa
sij=sij*1e-6;
% Plot results 
figure(1);
% Case (1)
sij1=sij(1,:);
plot(t,sij1,'k');
hold on;
% Case (2)
sij2=sij(2,:);
plot(t,sij2,'b');
% Case (3)
sij3=sij(3,:);
plot(t,sij3,'r');
% Case (4)
sij4=sij(4,:);
plot(t,sij4,'g');
% Case (5)
sij5=sij(5,:);
plot(t,sij5,'y');
% Case (6)
sij6=sij(6,:);
plot(t,sij6,'m');
hold off;