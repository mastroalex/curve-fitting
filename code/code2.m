clear all; close all; clc
%%
load mistery.mat


%% 
% ottimali in opportune variabili 
firstSignal=200;
lastSignal=40200;
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
a_temp=a(a>0);
shape_temp=shape(a>0);
clear a
clear shape
a=a_temp;
shape=shape_temp;
clear shape_temp
clear a_temp
diam=G*(a.*1e6).^(1/3); % electric parameters prop to a^(1/3) [um]

d_nominal=[5.2,6,7]; % nominal bead diameter [um]
n_bead=length(d_nominal);
diam_lim=[4.5 10.5];
shape_lim=[0.15 0.3];
vel_lim=[0.1 0.5];
n_bin=n_signal/10;
Color_orange='#D95319';
Color_blue='#4DBEEE';
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
% scatter plot
scatter_fig=figure();
scatter(diam,shape)
xlabel('Electric diameter [\mu m]')
ylabel('Shape parameters')
xlim(diam_lim)
ylim(shape_lim)
%%
% L=40e-6; %[mm] L=40 [um]
% velocity=L./(delta*1e-3); % where L [um] and delta [ms]
% velocity_fig=figure();
% scatter(diam,velocity)
%  % ATTENZIONE MANC QUALCOSA NELLE UNIT° DI MISURA !!!!
% xlabel('Electric diameter [\mu m]')
% ylabel('Velocity [mm/s]')
% xlim(diam_lim)
% ylim(vel_lim)
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

histogram_corrected2_fig=figure;
for i=1:n_bead
hold on
histogram(corrected_separate_value{i}(:,1),floor(n_bin/4),'EdgeColor',mycolor{i},'FaceColor',mycolor{i})
end
xlim(diam_lim)
xlabel('Corrected electric diameter [\mu m]')
ylabel('Count')
title('Corrected electric diameter')



%% extra test


n_bin_d=300;
n_bin_sig=300;
% limits from Errico_Caselli_SAB_2017 

electric_D_lim     =[4.5,10.5];
electric_D_norm_lim=[0.8,1.6];
velocity_lim       =[0.1,0.5];
sig_ov_del_lim     =[0.15,0.30];
figure()
histogram2(diam,shape,'DisplayStyle','tile','ShowEmptyBins','on', ...
    'NumBins',[n_bin_d,n_bin_sig],'XBinLimits',electric_D_lim,'YBinLimits',sig_ov_del_lim);
colorbar

figure()
histogram2(diam_corr,y_data,'DisplayStyle','tile','ShowEmptyBins','on', ...
    'NumBins',[n_bin_d,n_bin_sig],'XBinLimits',electric_D_lim,'YBinLimits',sig_ov_del_lim);
colorbar



%% Figures export
% insert path
% here is used path from Current Folder referencing path
path='figs2/';
exportgraphics(figure(signal_visualization_fig),strcat(path,'signal_visualization_fig','.pdf'),'BackgroundColor','none','ContentType','vector');
exportgraphics(figure(histogram_figs),strcat(path,'histogram_figs','.pdf'),'BackgroundColor','none','ContentType','vector');
exportgraphics(figure(scatter_fig),strcat(path,'scatter_fig','.pdf'),'BackgroundColor','none','ContentType','vector');
exportgraphics(figure(normalized_diam_figs),strcat(path,'normalized_diam_figs','.pdf'),'BackgroundColor','none','ContentType','vector');
exportgraphics(figure(scatter_corrected_fig),strcat(path,'scatter_corrected_fig','.pdf'),'BackgroundColor','none','ContentType','vector');
exportgraphics(figure(histogram_corrected_fig),strcat(path,'histogram_corrected_fig','.pdf'),'BackgroundColor','none','ContentType','vector');

exportgraphics(figure(scatter_corrected2_fig),strcat(path,'scatter_corrected2_fig','.pdf'),'BackgroundColor','none','ContentType','vector');
exportgraphics(figure(histogram_corrected2_fig),strcat(path,'histogram_corrected2_fig','.pdf'),'BackgroundColor','none','ContentType','vector');



