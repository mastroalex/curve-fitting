% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Mastrofini Alessandro
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Medical Engineering - University of Rome Tor Vergata
% Physiological Systems Modeling and Simulation
% F. Caselli, MSSF A.Y. 2021/2022
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Curve fitting for impedance micro cytometers
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

clear all; close all; clc
%%
load mistery.mat
%% Fitting the signals
% select range for signal selection
firstSignal=200;
lastSignal=35200;
n_signal=lastSignal-firstSignal;% total no. of signal considered
% initialize array to collect fitting parameters inside the loop for each
% signal considered
sigma=zeros(n_signal,1);
delta=sigma;
t_c=sigma;
a=sigma;
% fitting use a bipolar gaussian template like that:
% a*(exp(-((t-(t_c-delta/2)).^2/(2*sigma.^2)))-exp(-((t-(t_c+delta/2)).^2/(2*sigma.^2))))
fs=115e3; % sampling frequencies use for time scaling
tic
parfor i=1:(n_signal)
    j=i+firstSignal-1; % set the index to the corresponding signal
    % ectract fitting coefficients
    fitted=fitMySignal(mistery_data,j,'no');
    sigma(i)=fitted.sigma;
    delta(i)=fitted.delta;
    t_c(i)=fitted.t_c;
    a(i)=fitted.a;
    % scale back coefficient
    % fitting function normalize data to not work with small numbers
    % but use numbers around 1. So it is necessary to scale it back
    a(i)=a(i)*max(abs(mistery_data{j})); % scale amplitude with max
    Ns=length(mistery_data{j}); % number of sample
    time_step=1e3*(0:Ns-1)/fs;  % time step from [s] to [ms]
    % scale delta and t_c with last time % scale μs 
    delta(i)=delta(i)*time_step(end);
    t_c(i)=t_c(i)*time_step(end);
    sigma(i)=sigma(i)*time_step(end);
end
compute_time=toc; % save camputational cost
disp(['Completed in ', num2str(compute_time),' s'])

%% Plot collected data

% compute useful parameters

shape=sigma./delta; % shape parameters
% avoid numerical problems by selecting only positive max amplitude signal
a_temp=a(a>0);
shape_temp=shape(a>0);
delta_temp=delta(a>0);
clear a
clear shape
clear delta
a=a_temp;
shape=shape_temp;
delta=delta_temp;
clear shape_temp
clear a_temp
clear delta_temp
G=10.5; % From Errico [um / uA^(1/3)]
% electric diameters is proportional to a^(1/3) and is in [um]
% scale a into [uA] form previolsy scaled values
diam=G*(a.*1e6).^(1/3); % electric parameters prop to a^(1/3) [um]

% useful parameters
d_nominal=[5.2,6,7]; % nominal bead diameter [um]
n_bead=length(d_nominal);
% plotting limits
diam_lim=[4.5 10.5];
shape_lim=[0.15 0.3];
vel_lim=[0.1 0.5];
n_bin=n_signal/10; % histogram bins
% colours for plotting
Color_orange='#D95319';
Color_blue='#4DBEEE';
Color_green='#77AC30';
mycolor={Color_orange,Color_blue,Color_green};

% histogram 
histogram_figs=figure()
histogram(diam,n_bin,'EdgeAlpha',0.2);
title('Electrical diameter distribution') 
xlim(diam_lim)
ylabel('Count')
xlabel('Electrical diameter [\mu m]')

% scatter plot electric diameters vs shape parameters
scatter_fig=figure();
scatter(diam,shape)
xlabel('Electric diameter [\mu m]')
ylabel('Shape parameters')
xlim(diam_lim)
ylim(shape_lim)

% scatter plot electric diameters vs velocity
% velocity is defined as v=L/delta
% where L is the electrode interdistance
L=40e-6; %[mm] where L=40 [um]
% calculate and plot
velocity=L./(delta*1e-3); % where L [um] and delta [ms]
velocity_fig=figure();
scatter(diam,velocity)
xlabel('Electric diameter [\mu m]')
ylabel('Velocity [m/s]')
xlim(diam_lim)
ylim(vel_lim)

% due to the high numbers of singal could be useful to use density plot and not scatter ones

n_bin_d=300;
n_bin_sig=300;
% limits from Errico_Caselli_SAB_2017 
electric_D_lim     =[4.5,10.5];
electric_D_norm_lim=[0.8,1.6];
velocity_lim       =[0.1,0.5];
sig_ov_del_lim     =[0.15,0.30];

density_diam_shape=figure()
histogram2(diam,shape,'DisplayStyle','tile','ShowEmptyBins','on', ...
    'NumBins',[n_bin_d,n_bin_sig],'XBinLimits',electric_D_lim,'YBinLimits',sig_ov_del_lim);
colorbar
xlabel('Electric diameter [\mu m]')
ylabel('Shape parameters')

density_diam_velocity=figure()
histogram2(diam,velocity,'DisplayStyle','tile','ShowEmptyBins','on', ...
    'NumBins',[n_bin_d,n_bin_sig],'XBinLimits',electric_D_lim,'YBinLimits',velocity_lim);
colorbar
xlabel('Electric diameter [\mu m]')
ylabel('Velocity [m/s]')
%% Families selection and fitting

% remember to select from left to right for growing diameters of the beads
message = {'Click Ok and than select with the polygon the 3 different families from left'};
f = warndlg(message,'Warning');

% loop to select different families form the scatter plot
for i=1:n_bead
[X,Y]=getline(scatter_fig);
% manually select different indexes
[index_value]=inpolygon(diam,shape,X,Y);
% isolated the selected value for diameters and shape parameters
% with a char array of two column
selected_value{:,i}=[diam(index_value),shape(index_value)];
% compute the normalized diameters respect the relative nominal diameter
normalized_D{:,i}=diam(index_value)./d_nominal(i);
% compute fitting only for graphical purpose to give feedback in fig window
% the fitting is from diameter and shape parameters
fit_line(:,i)=polyfit(selected_value{i}(:,1),selected_value{i}(:,2),1)';
fitted_y=polyval(fit_line(:,i)',selected_value{i}(:,1));
fitted_value{i}=fitted_y;
% plot fitted line inside the scatter plot
hold on
plot(selected_value{i}(:,1),fitted_value{i},'LineWidth',2,'Color',mycolor{i})
% extract correspondig velocity values
velocity_corrected{i}=velocity(index_value);
clear X Y
end

% plot normalized value for the different families
normalized_diam_figs=figure;
hold on
for i=1:n_bead
plot(normalized_D{i},selected_value{i}(:,2),'*')
end

% computed fitting for compensation procedure
% fitting all the normalized diameters and shape parameters
% to extract only 2 coefficient for all the families
y_data=[];
x_data=[];
vel_data=[];
for i=1:n_bead
y_data=[y_data; selected_value{i}(:,2)];
x_data=[x_data; normalized_D{i}];
vel_data=[vel_data; velocity_corrected{i}];
end

fitting=polyfit(x_data,y_data,1);
fitting_valuated=polyval(fitting,x_data);
figure(normalized_diam_figs)
hold on
plot(x_data,fitting_valuated)

normalized_D=figure()
histogram2(x_data,y_data,'DisplayStyle','tile','ShowEmptyBins','on', ...
    'NumBins',[n_bin_d,n_bin_sig],'XBinLimits',[0.6, 1.6],'YBinLimits',sig_ov_del_lim);
colorbar
hold on
plot(x_data,fitting_valuated,'LineWidth',2)
xlabel('Normalized electric diameter [{\mu} m]')
ylabel('Shape parameters')
xlim([0.6, 1.6])
%% Compensation

% extract the coefficient of y=p1 x + p2 
p1=fitting(1);
p2=fitting(2);
% revert the coefficient to express normD=c1+c2.*shape
c1=-(p2/p1);
c2=1/p1;

% join diameter of selected value inside a single array
selected_diameter=[];
for i=1:n_bead
selected_diameter=[selected_diameter; selected_value{i}(:,1)];
end

% correct all the diameters with compensation formula
diam_corr=selected_diameter./(c1+c2*y_data);

density_diam_corr_shape=figure()
histogram2(diam_corr,y_data,'DisplayStyle','tile','ShowEmptyBins','on', ...
    'NumBins',[n_bin_d,n_bin_sig],'XBinLimits',electric_D_lim,'YBinLimits',sig_ov_del_lim);
colorbar
xlabel('Corrected electric diameter [\mu m]')
ylabel('Shape parameters')

density_diam_corr_velocity=figure()
histogram2(diam_corr,vel_data,'DisplayStyle','tile','ShowEmptyBins','on', ...
    'NumBins',[n_bin_d,n_bin_sig],'XBinLimits',electric_D_lim,'YBinLimits',velocity_lim);
colorbar
xlabel('Corrected electric diameter [\mu m]')
ylabel('Velocity [m/s]')

%% Plot corrected value separating different families

% compute separated corrected diametes
corrected_separate_value={};
for i=1:n_bead
corrected_separate_value{i}=[selected_value{i}(:,1)./(c1+c2*selected_value{i}(:,2)) , selected_value{i}(:,2)];
end

% histogram
histogram_corrected_fig=figure;
% fit the histrogram with a gaussian distribution
myhist={};
for i=1:n_bead
    hold on
    % single plot
    % histogram(corrected_separate_value{i}(:,1),floor(n_bin/4),'EdgeColor',mycolor{i},'FaceColor',mycolor{i})
    % the single plot is substituted from the fitting ones
    myhist{i}=histfit(corrected_separate_value{i}(:,1),floor(n_bin/4),'normal');
    % set color and transparency
    set(myhist{i}(1),'FaceColor',mycolor{i});
    set(myhist{i}(2),'color',mycolor{i})
    set(myhist{i}(2),'LineWidth',2)
    set(myhist{i}(2),'LineStyle','-')
    set(myhist{i}(1),'FaceAlpha',.4);
    % extract distribution parameters
    distribution_fit{i}=fitdist(corrected_separate_value{i}(:,1),'normal');
end
% graphics
xlim(diam_lim)
xlabel('Corrected electric diameter [\mu m]')
ylabel('Count')
title('Corrected electric diameter')
for i=1:n_bead
    j=i*2;
    myLegend{j-1}=strcat(string(d_nominal(i)),' {\mu}m beads');
    myLegend{j}=strcat('\sigma= '," ",string(distribution_fit{i}.sigma),', \mu= '," ",string(distribution_fit{i}.mu));
end
legend(myLegend)


%% Figures export
% insert path
% here is used path from Current Folder referencing path
path='figs2/';
exportgraphics(figure(signal_visualization_fig),strcat(path,'signal_visualization_fig','.pdf'),'BackgroundColor','none','ContentType','vector');
exportgraphics(figure(histogram_figs),strcat(path,'histogram_figs','.pdf'),'BackgroundColor','none','ContentType','vector');
exportgraphics(figure(scatter_fig),strcat(path,'scatter_fig','.pdf'),'BackgroundColor','none','ContentType','vector');
exportgraphics(figure(histogram_figs),strcat(path,'histogram_data','.pdf'),'BackgroundColor','none','ContentType','vector');

exportgraphics(figure(density_diam_shape),strcat(path,'density_diam_shape','.pdf'),'BackgroundColor','none','ContentType','vector');
exportgraphics(figure(density_diam_velocity),strcat(path,'density_diam_velocity','.pdf'),'BackgroundColor','none','ContentType','vector');
exportgraphics(figure(normalized_D),strcat(path,'normalized_D','.pdf'),'BackgroundColor','none','ContentType','vector');


exportgraphics(figure(density_diam_corr_velocity),strcat(path,'density_diam_corr_velocity','.pdf'),'BackgroundColor','none','ContentType','vector');
exportgraphics(figure(density_diam_corr_shape),strcat(path,'density_diam_corr_shape','.pdf'),'BackgroundColor','none','ContentType','vector');

exportgraphics(figure(histogram_corrected_fig),strcat(path,'histogram_corrected_fig','.pdf'),'BackgroundColor','none','ContentType','vector');



