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

%% Signal visualization
% select index of the signal to plot
% from 1 to 57361
signal_to_plot=[5 10 400 700 35698];
% different marks for different signal
all_marks = {'o','+','*','x','^','v','>','p','h'};
% setup figure and plot selected signal
signal_visualization_fig=figure();
hold on
for i=1:length(signal_to_plot)
    % update the marker from the list
    marker=all_marks{mod(i,length(all_marks))};
    % plot the signal of index selected from signal_to_plot
    % the amplitude is scaled with 1e6 to scale from [A] to [uA]
    plot(1e6*mistery_data{signal_to_plot(i)},marker)
end
ylabel('I_{diff} [\mu A]')
xlabel('Sample')
grid on
%% Fitting test
% use like fitMySignal(index,mistery_data,WantPlot)
% WantPlot to 'no' and no plot
% WantPlot to 'yes' return also plot figure
% default is set to 'yes' to obtaine plotted figure
% the function return plot figure and fitting output
index=100;

%fitted=fitMySignal(mistery_data,index,'no');
fitted=fitMySignal(mistery_data,index);

%% Fitting complete
% select range for signal selection
firstSignal=400;
lastSignal=800;
%lastSignal=57361;
n_signal=lastSignal-firstSignal; % total no. of signal considered
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
    % scale delta and t_c with last time
    delta(i)=delta(i)*time_step(end);
    t_c(i)=t_c(i)*time_step(end);
    sigma(i)=sigma(i)*time_step(end);
end
compute_time=toc; % save camputational cost
disp(['Completed in ', num2str(compute_time),' s'])

%% Plot collected data

% compute useful parameters
shape=sigma./delta; % shape parameters
%G=1;  % electric gain [um / uA^(1/3)]
G=10.5; % From Errico et al [um / uA^(1/3)]
% electric diameters is proportional to a^(1/3) and is in [um]
% scale a into [uA] form previolsy scaled values
diam=G*(a.*1e6).^(1/3);

% useful parameters for plotting
d_nominal=[5.2,6,7]; % nominal bead diameter [um]
n_bead=length(d_nominal); % total n of beads
% plotting limits
diam_lim=[4.5 10.5];
shape_lim=[0.15 0.3];
vel_lim=[0.1 0.5];
% histogram bins
n_bin=n_signal/4;
% colours for plotting
Color_orange='#D95319';
Color_blue='#0072BD';
Color_green='#77AC30';
mycolor={Color_orange,Color_blue,Color_green};

% histogram
histogram_figs=figure()
histogram(diam,n_bin);
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
%% Families selection and fitting

% remember to select from left to right for growing diameters of the beads
message = {'Click Ok and than select with the polygon the 3 different families from left'};
f = warndlg(message,'Warning');

% loop to select different families
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
    % the fitting is for diameter and shape parameters
    fit_line(:,i)=polyfit(selected_value{i}(:,1),selected_value{i}(:,2),1)';
    fitted_y=polyval(fit_line(:,i)',selected_value{i}(:,1));
    fitted_value{i}=fitted_y;
    hold on
    % plot fitted line inside the scatter plot
    plot(selected_value{i}(:,1),fitted_value{i},'LineWidth',2,'Color',mycolor{i})
    clear X Y
    % extract correspondig velocity values
    velocity_corrected{i}=velocity(index_value);
end

% plot normalized value for the different families
normalized_diam_figs=figure;
hold on
for i=1:n_bead
    plot(normalized_D{i},selected_value{i}(:,2),'*','MarkerEdgeColor',mycolor{i})
end

% computed fitting for compensation procedure
% fitting all the normalized diameters and shape parameters
% to extract only 2 coefficient for all the families
y_data=[];
x_data=[];
for i=1:n_bead
    % append data inside a single array
    y_data=[y_data; selected_value{i}(:,2)];
    x_data=[x_data; normalized_D{i}];
end
% it is possibile to use polyfit with linear fitting
% % plot resulting value
% %fitting=polyfit(x_data,y_data,1);
% %fitting_valuated=polyval(fitting,x_data);
% % figure(normalized_diam_figs)
% % hold on
% % plot(x_data,fitting_valuated)

% but is is possible also to use linear regression model
% and so extract directly R2 and coefficients
fitting=fitlm(x_data,y_data);
figure(normalized_diam_figs)
hold on
% extract the coefficient of y=p2 x + p1 
p1=table2array(fitting.Coefficients(1,1));
p2=table2array(fitting.Coefficients(2,1));
% plot value 
plot(x_data,p2.*x_data+p1,'-k')
legend([strcat(string(d_nominal),' {\mu}m beads'), strcat('r^2 = ',num2str(fitting.Rsquared.Ordinary,2))])
xlabel('Normalized electric diameters')
ylabel('Shape parameter')
%% Compensation procedure
% extract the coefficient of y=p2 x + p1 
p1=table2array(fitting.Coefficients(1,1));
p2=table2array(fitting.Coefficients(2,1));
% revert the coefficient to express normD=c1+c2.*shape
c1=-(p1/p2);
c2=1/p2;

% join diameter of selected value inside a single array
selected_diameter=[];
for i=1:n_bead
    selected_diameter=[selected_diameter; selected_value{i}(:,1)];
end

% correct all the diameters with compensation formula
diam_corr=selected_diameter./(c1+c2*y_data);

% plot corrected value inside a single families
% scatter plot
scatter_corrected_fig=figure;
scatter(diam_corr,y_data)
xlabel('Corrected electric diameter [\mu m]')
ylabel('Shape parameters')
xlim(diam_lim)
ylim(shape_lim)
% histrogram t
histogram_corrected_fig=figure;
histogram(diam_corr,n_bin)
xlabel('Corrected electric diameter [\mu m]')
ylabel('Count')
title('Corrected electric diameter')
xlim(diam_lim)

%% Plot corrected value separating different families

% compute separated corrected diametes
corrected_separate_value={};
for i=1:n_bead
    corrected_separate_value{i}=[selected_value{i}(:,1)./(c1+c2*selected_value{i}(:,2)) , selected_value{i}(:,2)];
end

% plot scatter 
scatter_corrected2_fig=figure;
for i=1:n_bead
    hold on
    scatter(corrected_separate_value{i}(:,1),corrected_separate_value{i}(:,2),'MarkerEdgeColor',mycolor{i})
end
xlabel('Corrected electric diameter [\mu m]')
ylabel('Shape parameters')
xlim(diam_lim)
ylim(shape_lim)

% histogram
histogram_corrected2_fig=figure;
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
    set(myhist{i}(1),'FaceAlpha',.2);
    % extract distribution parameters
    distribution_fit{i}=fitdist(corrected_separate_value{i}(:,1),'normal');
end
% plot graphics
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

% velocity scatter plot
velocity_corrected_fig=figure();
for i=1:n_bead
    hold on
    scatter(corrected_separate_value{i}(:,1),velocity_corrected{i},'MarkerEdgeColor',mycolor{i})
end
xlabel('Corrected electric diameter [\mu m]')
ylabel('Velocity [m/s]')
xlim(diam_lim)
ylim(vel_lim)

%% Figures export
% insert path
% here is used path from Current Folder referencing path
path='figs/';
%exportgraphics(figure(signal_visualization_fig),strcat(path,'signal_visualization_fig','.pdf'),'BackgroundColor','none','ContentType','vector');
exportgraphics(figure(histogram_figs),strcat(path,'histogram_figs','.pdf'),'BackgroundColor','none','ContentType','vector');
exportgraphics(figure(scatter_fig),strcat(path,'scatter_fig','.pdf'),'BackgroundColor','none','ContentType','vector');
exportgraphics(figure(velocity_fig),strcat(path,'velocity_fig','.pdf'),'BackgroundColor','none','ContentType','vector');
exportgraphics(figure(normalized_diam_figs),strcat(path,'normalized_diam_figs','.pdf'),'BackgroundColor','none','ContentType','vector');
exportgraphics(figure(scatter_corrected_fig),strcat(path,'scatter_corrected_fig','.pdf'),'BackgroundColor','none','ContentType','vector');
exportgraphics(figure(histogram_corrected_fig),strcat(path,'histogram_corrected_fig','.pdf'),'BackgroundColor','none','ContentType','vector');

exportgraphics(figure(scatter_corrected2_fig),strcat(path,'scatter_corrected2_fig','.pdf'),'BackgroundColor','none','ContentType','vector');
exportgraphics(figure(histogram_corrected2_fig),strcat(path,'histogram_corrected2_fig','.pdf'),'BackgroundColor','none','ContentType','vector');
exportgraphics(figure(velocity_corrected_fig),strcat(path,'velocity_corrected_fig','.pdf'),'BackgroundColor','none','ContentType','vector');

save('final_test_workspace.mat')