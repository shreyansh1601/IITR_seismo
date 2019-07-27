% Loading from file and visualizing density and enthalpy maps
% computed for pirolite mantle (files m895_ro, m895_hh) and
% MORB oceanic crust (files morn_ro, morn_hh)
% Function Melt_fraction()
% This function compute melt fraction (xmelt) and latent heat (hlat)
% at given pressure (ppa), temperature (mtk)and rock type (rock)
% Function returns solution for new temperature (vx,vy,pr)
% and distribution of residuals (resx,resy,resc)
function[tknum,pbnum,tkmin,pbmin,tkstp,pbstp,td]=loading_database()


% Loading density data for the mantle
% Opening data file
fdata=fopen('m895_ro','rt');
% Skipping rock title
a=fscanf(fdata,'%s',1);
% Loading resolutions for T and P
tknum=fscanf(fdata,'%f',1);
pbnum=fscanf(fdata,'%f',1);
% Loading starting values T(K) and P(bar)
tkmin=fscanf(fdata,'%f',1);
pbmin=fscanf(fdata,'%f',1);
% Loading steps for T(K) and P(bar)
tkstp=fscanf(fdata,'%f',1);
pbstp=fscanf(fdata,'%f',1);
% Skipping titles for P and T
b=fscanf(fdata,'%s',2);
% Creating database array
td=zeros(2,2,pbnum,tknum);
% Loading density data
for i=1:1:pbnum
    for j=1:1:tknum
        td(1,1,i,j)=fscanf(fdata,'%f',1);
    end
end
% Closing data file
fclose(fdata);

% Loading density data for the mantle
% Opening data file
fdata=fopen('m895_hh','rt');
% Skipping rock title
a=fscanf(fdata,'%s',1);
% Loading resolutions for T and P
tknum=fscanf(fdata,'%f',1);
pbnum=fscanf(fdata,'%f',1);
% Loading starting values T(K) and P(bar)
tkmin=fscanf(fdata,'%f',1);
pbmin=fscanf(fdata,'%f',1);
% Loading steps for T(K) and P(bar)
tkstp=fscanf(fdata,'%f',1);
pbstp=fscanf(fdata,'%f',1);
% Skipping titles for P and T
b=fscanf(fdata,'%s',2);
% Loading density data
for i=1:1:pbnum
    for j=1:1:tknum
        td(1,2,i,j)=fscanf(fdata,'%f',1);
    end
end
% Closing data file
fclose(fdata);


% Loading density data for the crust
% Opening data file
fdata=fopen('morn_ro','rt');
% Skipping rock title
a=fscanf(fdata,'%s',1);
% Loading resolutions for T and P
tknum=fscanf(fdata,'%f',1);
pbnum=fscanf(fdata,'%f',1);
% Loading starting values T(K) and P(bar)
tkmin=fscanf(fdata,'%f',1);
pbmin=fscanf(fdata,'%f',1);
% Loading steps for T(K) and P(bar)
tkstp=fscanf(fdata,'%f',1);
pbstp=fscanf(fdata,'%f',1);
% Skipping titles for P and T
b=fscanf(fdata,'%s',2);
% Loading density data
for i=1:1:pbnum
    for j=1:1:tknum
        td(2,1,i,j)=fscanf(fdata,'%f',1);
    end
end
% Closing data file
fclose(fdata);

% Loading density data for the mantle
% Opening data file
fdata=fopen('morn_hh','rt');
% Skipping rock title
a=fscanf(fdata,'%s',1);
% Loading resolutions for T and P
tknum=fscanf(fdata,'%f',1);
pbnum=fscanf(fdata,'%f',1);
% Loading starting values T(K) and P(bar)
tkmin=fscanf(fdata,'%f',1);
pbmin=fscanf(fdata,'%f',1);
% Loading steps for T(K) and P(bar)
tkstp=fscanf(fdata,'%f',1);
pbstp=fscanf(fdata,'%f',1);
% Skipping titles for P and T
b=fscanf(fdata,'%s',2);
% Loading density data
for i=1:1:pbnum
    for j=1:1:tknum
        td(2,2,i,j)=fscanf(fdata,'%f',1);
    end
end
% Closing data file
fclose(fdata);


% Defining new figure
figure(4);

% 2D colormap
% Creating a vector of arguments for two axes
for i=1:1:pbnum
    p(i)=pbmin*1e-4+(i-1)*pbstp*1e-4; % Pressure GPa (1 bar = 10^5 Pa = 10^(-4) GPa
end
for j=1:1:tknum
    t(j)=tkmin-273+(j-1)*tkstp; % T C
end

% Plotting mantle density
c=zeros(pbnum,tknum);
for i=1:1:pbnum
    for j=1:1:tknum
        c(i,j)=td(1,1,i,j);
    end
end
subplot(2,2,1);     % defining 2nd plotting area in the 3x1 figure
pcolor(t,p,c);      % making colormap
box on;             % making box around the plot
title('Pirolite density, kg/m^3');   % title of the plot
xlabel('T,C');        % title of horizontal axis
ylabel('P,GPa');        % title of vertical axis
shading interp;     % making smooth transitions between colors
colorbar;           % showing colorbar for the map

% Plotting crust density
c=zeros(pbnum,tknum);
for i=1:1:pbnum
    for j=1:1:tknum
        c(i,j)=td(2,1,i,j);
    end
end
subplot(2,2,2);     % defining 2nd plotting area in the 3x1 figure
pcolor(t,p,c);      % making colormap
box on;             % making box around the plot
title('MORB density, kg/m^3');   % title of the plot
xlabel('T,C');        % title of horizontal axis
ylabel('P,GPa');        % title of vertical axis
shading interp;     % making smooth transitions between colors
colorbar;           % showing colorbar for the map

% Plotting mantle enthalpy
c=zeros(pbnum,tknum);
for i=1:1:pbnum
    for j=1:1:tknum
        c(i,j)=td(1,2,i,j);
    end
end
subplot(2,2,3);     % defining 2nd plotting area in the 3x1 figure
pcolor(t,p,c);      % making colormap
box on;             % making box around the plot
title('Pirolite enthalpy, J/kg');   % title of the plot
xlabel('T,C');        % title of horizontal axis
ylabel('P,GPa');        % title of vertical axis
shading interp;     % making smooth transitions between colors
colorbar;           % showing colorbar for the map

% Plotting crust enthalpy
c=zeros(pbnum,tknum);
for i=1:1:pbnum
    for j=1:1:tknum
        c(i,j)=td(2,2,i,j);
    end
end
subplot(2,2,4);     % defining 2nd plotting area in the 3x1 figure
pcolor(t,p,c);      % making colormap
box on;             % making box around the plot
title('MORB enthalpy, J/kg');   % title of the plot
xlabel('T,C');        % title of horizontal axis
ylabel('P,GPa');        % title of vertical axis
shading interp;     % making smooth transitions between colors
colorbar;           % showing colorbar for the map


