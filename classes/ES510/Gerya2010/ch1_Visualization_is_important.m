% visualization of ‘sin’ and ‘cos’ functions with ‘plot’, ‘pcolor’, ‘contour’ and  ‘surf’

% Clearing all variables and arrays
clear all;
% Clearing figures
clf;

% Defining new figure
figure(1)

% 2D plot
% Creating a vector of arguments for horizontal axis
x=0:1:100;
% Creating a function to be plotted for vertical axis
y=sin(x/5).*cos(x/10);
% Plotting the function
subplot(1,3,1);     % defining 1st plotting area in the 3x1 figure
plot(x,y,'k . -.'); % making plot: 'k' = black color, '.' = dotted markers, '-.' = dotted-dashed line
box on;             % making box around the plot
title('2D plot');   % title of the plot
xlabel('x');        % title of horizontal axis
ylabel('y=sin(x/5)*cos(x/10)'); % title of vertical axis



% 2D colormap
% Creating a vector of arguments for two axes
x=0:1:100;
y=0:1:100;
% Creating a function to be plotted for vertical axis
% Smallest index in an array is 1 (not 0)
for i=1:1:101
    for j=1:1:101
        z(i,j)=sin(x(j)/20).*cos(y(i)/10);
    end
end
% Plotting the function
subplot(1,3,2);     % defining 2nd plotting area in the 3x1 figure
pcolor(x,y,z);      % making colormap
box on;             % making box around the plot
title('2D colormap');   % title of the plot
xlabel('x');        % title of horizontal axis
ylabel('y');        % title of vertical axis
shading interp;     % making smooth transitions between colors
colorbar;           % showing colorbar for the map
hold on;            % continuing plotting on the colormap 
[C, h] = contour (x,y,z,'k');   % drawing isolines
clabel(C,h,'Color','k');        % labelling isolines
hold off;           % stop plotting on the colormap



% 3D surface
% Plotting the previous colormap function
subplot(1,3,3);     % defining 3nd plotting area in the 3x1 figure
surf(x,y,z);        % making color surface
title('3D surface');% title of the plot
xlabel('x');        % title of horizontal axis
ylabel('y');        % title of lateral axis
zlabel('z=sin(x/20)*cos(y/10)'); % title of vertical axis
shading interp;     % making smooth transitions between colors
colorbar;           % showing colorbar for the colored surface
light;              % shedding light to the surface
lighting phong;     % making nicely reflecting surface


