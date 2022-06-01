clear all; close all; clc
%%
load mistery.mat

%% Signal visualization
% select index of the signal to plot
% from 1 to 57361
signal_to_plot=[5 10 400 700 35698];
% different marks for different signal
all_marks = {'o','+','*','x','^','v','>','p','h'};
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
%% Fitting (from previous lesson)
% use like fitMySignal(index,mistery_data,WantPlot)
% WantPlot to 'no' and no plot
% WantPlot to 'yes' return also plot figure
% default is set to 'yes' to obtaine plotted figure
% the function return plot figure and fitting output
index=100;

%fitted=fitMySignal(mistery_data,index,'no');
fitted=fitMySignal(mistery_data,index);

%% 
% Predisporre un ciclo for su un tot di segnali, memorizzando i parametri
% ottimali in opportune variabili 
firstSignal=400;
lastSignal=800;
%lastSignal=57361;
n_signal=lastSignal-firstSignal;
sigma=zeros(lastSignal-firstSignal,1);
delta=sigma;
t_c=sigma;
a=sigma;
fs=115e3;
tic
parfor i=1:(lastSignal-firstSignal)
    j=i+firstSignal-1;
    fitted=fitMySignal(mistery_data,j,'no');
    sigma(i)=fitted.sigma;
    delta(i)=fitted.delta;
    t_c(i)=fitted.t_c;
    a(i)=fitted.a;
    % scale back
    a(i)=a(i)*max(abs(mistery_data{j})); % scale amplitude with max
    Ns=length(mistery_data{j});
    time_step=1e3*(0:Ns-1)/fs;  % scale delta and t_c with last time % scale μs 
    delta(i)=delta(i)*time_step(end);
    t_c(i)=t_c(i)*time_step(end);
    sigma(i)=sigma(i)*time_step(end);
end
compute_time=toc;
disp(['Completed in ', num2str(compute_time),' s'])

%% PLOT
shape=sigma./delta; % shape parameters
%G=1;                % electric gain [um / uA^(1/3)]
G=10.5; % From Errico [um / uA^(1/3)]
% a is [A] but I need it in uA --> *1e6
diam=G*(a.*1e6).^(1/3); % electric parameters prop to a^(1/3) [um]
d_nominal=[5.2,6,7]; % nominal bead diameter [um]
n_bead=length(d_nominal);
diam_lim=[4.5 10.5];
shape_lim=[0.15 0.3];
vel_lim=[0.1 0.5];
n_bin=n_signal/4;
Color_orange='#D95319';
Color_blue='#0072BD';
Color_green='#77AC30';
mycolor={Color_orange,Color_blue,Color_green};
% histogram 

histogram_figs=figure()
%histogram(diam);
% to use M BIN
histogram(diam,n_bin);
title('Electrical diameter distribution') 
xlim(diam_lim)
ylabel('Count')
xlabel('Electrical diameter [\mu m]')
%% SCATTER PLOT and LINEAR FIT
% scatter plot
scatter_fig=figure();
scatter(diam,shape)
xlabel('Electric diameter [\mu m]')
ylabel('Shape parameters')
xlim(diam_lim)
ylim(shape_lim)
%%
L=40e-6; %[mm] L=40 [um]
velocity=L./(delta*1e-3); % where L [um] and delta [ms]
velocity_fig=figure();
scatter(diam,velocity)
 % ATTENZIONE MANC QUALCOSA NELLE UNIT° DI MISURA !!!!
xlabel('Electric diameter [\mu m]')
ylabel('Velocity [mm/s]')
xlim(diam_lim)
ylim(vel_lim)
%% 


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
% the fitting is from diameter and shape parameters
fit_line(:,i)=polyfit(selected_value{i}(:,1),selected_value{i}(:,2),1)';
fitted_y=polyval(fit_line(:,i)',selected_value{i}(:,1));
fitted_value{i}=fitted_y;
hold on
plot(selected_value{i}(:,1),fitted_value{i},'LineWidth',2,'Color',mycolor{i})
clear X Y
velocity_corrected{i}=velocity(index_value);
end

normalized_diam_figs=figure;
hold on
% plot normalized value for the different families
for i=1:n_bead
plot(normalized_D{i},selected_value{i}(:,2),'*')
end

% computed fitting for compensation procedure
y_data=[];
x_data=[];
for i=1:n_bead
y_data=[y_data; selected_value{i}(:,2)];
x_data=[x_data; normalized_D{i}];
end

fitting=polyfit(x_data,y_data,1);
fitting_valuated=polyval(fitting,x_data);
figure(normalized_diam_figs)
hold on
plot(x_data,fitting_valuated)

%% 
p1=fitting(1);
p2=fitting(2);
c1=-(p2/p1);
c2=1/p1;

selected_diameter=[];
for i=1:n_bead
selected_diameter=[selected_diameter; selected_value{i}(:,1)];
end

diam_corr=selected_diameter./(c1+c2*y_data);

scatter_corrected_fig=figure;
scatter(diam_corr,y_data)
xlabel('Corrected electric diameter [\mu m]')
ylabel('Shape parameters')
xlim(diam_lim)
ylim(shape_lim)

histogram_corrected_fig=figure;
histogram(diam_corr,n_bin)
xlabel('Corrected electric diameter [\mu m]')
ylabel('Count')
title('Corrected electric diameter')
xlim(diam_lim)

%% separate corrected value
corrected_separate_value={};
for i=1:n_bead
corrected_separate_value{i}=[selected_value{i}(:,1)./(c1+c2*selected_value{i}(:,2)) , selected_value{i}(:,2)];
end

scatter_corrected2_fig=figure;

for i=1:n_bead
hold on
scatter(corrected_separate_value{i}(:,1),corrected_separate_value{i}(:,2),'MarkerEdgeColor',mycolor{i})
end
xlabel('Corrected electric diameter [\mu m]')
ylabel('Shape parameters')
xlim(diam_lim)
ylim(shape_lim)

%%


histogram_corrected2_fig=figure;
myhist={};
for i=1:n_bead
hold on
%histogram(corrected_separate_value{i}(:,1),floor(n_bin/4),'EdgeColor',mycolor{i},'FaceColor',mycolor{i})
myhist{i}=histfit(corrected_separate_value{i}(:,1),floor(n_bin/4),'normal');
set(myhist{i}(1),'FaceColor',mycolor{i}); 
set(myhist{i}(2),'color',mycolor{i})
set(myhist{i}(2),'LineWidth',2)
set(myhist{i}(2),'LineStyle','-')
set(myhist{i}(1),'FaceAlpha',.2);
distribution_fit{i}=fitdist(corrected_separate_value{i}(:,1),'normal');
end
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

%%

velocity_corrected_fig=figure();

for i=1:n_bead
hold on
scatter(corrected_separate_value{i}(:,1),velocity_corrected{i},'MarkerEdgeColor',mycolor{i})
end
xlabel('Corrected electric diameter [\mu m]')
ylabel('Velocity [mm/s]')
xlim(diam_lim)
ylim(vel_lim)

%% Figures export
% insert path
% here is used path from Current Folder referencing path
path='figs/';
exportgraphics(figure(signal_visualization_fig),strcat(path,'signal_visualization_fig','.pdf'),'BackgroundColor','none','ContentType','vector');
exportgraphics(figure(histogram_figs),strcat(path,'histogram_figs','.pdf'),'BackgroundColor','none','ContentType','vector');
exportgraphics(figure(scatter_fig),strcat(path,'scatter_fig','.pdf'),'BackgroundColor','none','ContentType','vector');
exportgraphics(figure(normalized_diam_figs),strcat(path,'normalized_diam_figs','.pdf'),'BackgroundColor','none','ContentType','vector');
exportgraphics(figure(scatter_corrected_fig),strcat(path,'scatter_corrected_fig','.pdf'),'BackgroundColor','none','ContentType','vector');
exportgraphics(figure(histogram_corrected_fig),strcat(path,'histogram_corrected_fig','.pdf'),'BackgroundColor','none','ContentType','vector');

exportgraphics(figure(scatter_corrected2_fig),strcat(path,'scatter_corrected2_fig','.pdf'),'BackgroundColor','none','ContentType','vector');
exportgraphics(figure(histogram_corrected2_fig),strcat(path,'histogram_corrected2_fig','.pdf'),'BackgroundColor','none','ContentType','vector');
exportgraphics(figure(velocity_corrected_fig),strcat(path,'velocity_corrected_fig','.pdf'),'BackgroundColor','none','ContentType','vector');



