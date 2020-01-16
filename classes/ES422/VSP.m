% Problem Set1 : Question 3
% Vipul Silwal

clear all
close all

N = 100;
%N = 4;

z_inital = 0;
z_final = 20;

% discretizing the depth space
z_points = linspace(z_inital,z_final,N+1);
sensor_depth = z_points(2:length(z_points));

% creating mid-points for each layer
for i=1:length(z_points)-1
    z_mid(i) = z_points(i)+(z_points(i+1) - z_points(i))/2;
end

% making design matrix
G_base = tril(ones(N,N));
G = XXX

% making model_true space
for i=1:length(z_points)-1
    m_true(i,1) = XXX
end

% calculating data_true
d_true = G*m_true;

% evaluating data_exact
for i=1:length(z_points)-1
    temp = 0.2/((8*i)+996);
    if i>1
        d_exact(i) = temp + d_exact(i-1);
    else
        d_exact(i) = temp;
    end
end
d_exact = d_exact';

% plotting
figure;
subplot(2,1,1);
hold on
plot(z_points(2:length(z_points)),d_true,'xb');
plot(z_points(2:length(z_points)),d_exact,'or');
legend('d\_true','d\_exact');
xlabel('z (sensor depth)');
ylabel('data\_space (travel time)');
title('Data sapce');
grid on;

%======================Inversion=====================

% solving for m
m = XXX;

% plotting
subplot(2,1,2);
hold on
plot(z_mid,m_true,'xb');
plot(z_mid,m,'or');
legend('m\_true','m\_inv');
ylabel('model\_parameter (slowness)');
xlabel('z (interval midpoint)');
title('Model sapce');
grid on;

%===================Adding noise=====================

std = .05;
error = (std)*randn(N,1);
d_obs = d_exact+error;
m_inv = XXX;
figure

% plotting
subplot(2,1,1);
hold on
plot(z_points(2:length(z_points)),d_true,'xb');
plot(z_points(2:length(z_points)),d_obs,'or');
legend('d\_true','d\_exact+noise');
xlabel('z (sensor depth)');
ylabel('data\_space (travel time)');
title('Data sapce + noise');
grid on;

subplot(2,1,2);
hold on
plot(z_mid,m_true,'xb');
plot(z_mid,m_inv,'or');
legend('m\_true','m\_inv');
ylabel('model\_parameter (slowness)');
xlabel('z (interval midpoint)');
title('Inversion after adding noise to data');
grid on;

%============additional figures=====================
if (0)
    x=-20:.01:20;
    
    plot(x,[0],'-');
    hold on
    z = -20:.01:0;
    plot([0],z,'-');
    
    for i = 1:length(z)/100
        r=z(100*i);
        y=sqrt(r^2-x.^2);
        Y=real(y);
        plot(x,-Y);
    end
end
