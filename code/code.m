clear all; close all; clc
%%
load mistery.mat

%% Signal visualization
close all
figure()
hold on
starting_signal=1;
ending_signal=5;
for i=starting_signal:ending_signal
    plot(mistery_data{i},'*')
end
ylabel('I_{Diff} [A]')
xlabel('Distance along channel [\mu M]')
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
    time_step=1e3*(0:Ns-1)/fs;  % scale delta and t_c with last time 
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

diam_lim=[4.5 10.5];
shape_lim=[0.15 0.3];
vel_lim=[0.1 0.5];
n_bin=n_signal/4;
Color_orange='#D95319';
Color_blue='#0072BD';
Color_green='#77AC30';
% histogram 

figure()
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
figure();
scatter(diam,velocity)
 % ATTENZIONE MANC QUALCOSA NELLE UNIT° DI MISURA !!!!
xlabel('Electric diameter [\mu m]')
ylabel('Velocity [mm/s]')
xlim(diam_lim)
ylim(vel_lim)
%% 

message = {'Click Ok and than select with the polygon the 3 different families'};
f = warndlg(message,'Warning');

[X,Y]=getline(scatter_fig);
[index_value_family1]=inpolygon(diam,shape,X,Y);
selected_value_family1=[diam(index_value_family1),shape(index_value_family1)];
fit_line_family1=polyfit(selected_value_family1(:,1),selected_value_family1(:,2),1);
fitted_y_family1=polyval(fit_line_family1,selected_value_family1(:,1));
hold on
plot(selected_value_family1(:,1),fitted_y_family1,'LineWidth',2,'Color',Color_orange)
clear X Y
[X,Y]=getline(scatter_fig);
[index_value_family2]=inpolygon(diam,shape,X,Y);
selected_value_family2=[diam(index_value_family2),shape(index_value_family2)];
fit_line_family2=polyfit(selected_value_family2(:,1),selected_value_family2(:,2),1);
fitted_y_family2=polyval(fit_line_family2,selected_value_family2(:,1));
hold on
plot(selected_value_family2(:,1),fitted_y_family2,'LineWidth',2,'Color',Color_blue)

clear X Y
[X,Y]=getline(scatter_fig);
[index_value_family3]=inpolygon(diam,shape,X,Y);
selected_value_family3=[diam(index_value_family3),shape(index_value_family3)];
fit_line_family3=polyfit(selected_value_family3(:,1),selected_value_family3(:,2),1);
fitted_y_family3=polyval(fit_line_family3,selected_value_family3(:,1));
hold on
plot(selected_value_family3(:,1),fitted_y_family3,'LineWidth',2,'Color',Color_green)


%% 
% retta y=p1x+p2 -> x=y/p1 - p2/p1 -> - p2/p1 + y/p1
% c1+c2(shape)
p1_family1=fit_line_family1(1);
p2_family1=fit_line_family1(2);
c1_family1=-(p2_family1/p1_family1);
c2_family1=1/p1_family1;

diam_corr_family1=diam(index_value_family1)./(c1_family1+c2_family1*shape(index_value_family1));

p1_family2=fit_line_family2(1);
p2_family2=fit_line_family2(2);
c1_family2=-(p2_family2/p1_family2);
c2_family2=1/p1_family2;

diam_corr_family2=diam(index_value_family2)./(c1_family2+c2_family2*shape(index_value_family2));

p1_family3=fit_line_family3(1);
p2_family3=fit_line_family3(2);
c1_family3=-(p2_family3/p1_family3);
c2_family3=1/p1_family3;

diam_corr_family3=diam(index_value_family3)./(c1_family3+c2_family3*shape(index_value_family3));

figure()
histogram(diam_corr_family1,n_bin,'FaceColor',Color_orange);
hold on
histogram(diam_corr_family2,n_bin,'FaceColor',Color_blue);
histogram(diam_corr_family3,n_bin,'FaceColor',Color_green);
xlim(diam_lim)

figure()
scatter(diam_corr_family1,shape(index_value_family1),'Color',Color_orange)
xlim([0,2])

%% WHY THIS VALIUES ARE NORMALIZED ???
%% Test by numltiply for the mean
figure()
histogram(diam_corr_family1*mean(diam(index_value_family1)),n_bin,'FaceColor',Color_orange,'EdgeColor','none');
hold on
histogram(diam_corr_family2*mean(diam(index_value_family2)),n_bin,'FaceColor',Color_blue,'EdgeColor','none');
histogram(diam_corr_family3*mean(diam(index_value_family3)),n_bin,'FaceColor',Color_green,'EdgeColor','none');
xlim(diam_lim)

figure()
scatter(diam_corr_family1*mean(diam(index_value_family1)),shape(index_value_family1),'MarkerEdgeColor',Color_orange)
hold on
scatter(diam_corr_family2*mean(diam(index_value_family2)),shape(index_value_family2),'MarkerEdgeColor',Color_blue)
scatter(diam_corr_family3*mean(diam(index_value_family3)),shape(index_value_family3),'MarkerEdgeColor',Color_green)

xlim(diam_lim)




%% 
n_bin_d=200;
n_bin_sig=200;
% limits from Errico_Caselli_SAB_2017 

electric_D_lim     =[4.5,10.5];
electric_D_norm_lim=[0.8,1.6];
velocity_lim       =[0.1,0.5];
sig_ov_del_lim     =[0.15,0.30];
figure()
histogram2(diam_corr_family1,shape,'DisplayStyle','tile','ShowEmptyBins','on', ...
    'NumBins',[n_bin_d,n_bin_sig],'XBinLimits',electric_D_norm_lim,'YBinLimits',sig_ov_del_lim);
colorbar