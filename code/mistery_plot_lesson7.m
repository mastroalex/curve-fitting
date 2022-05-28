%% Lezione 7 - 22/03/2022
%%
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

%% Start lesson 7

%% STEP 0 
% Far girare il curve fitting su un segnale

% DONE

%% STEP 1 
% Predisporre un ciclo for su un tot di segnali, memorizzando i parametri
% ottimali in opportune variabili 
firstSignal=200;
lastSignal=400;
%lastSignal=57361;
sigma=zeros(lastSignal,1);
delta=sigma;
t_c=sigma;
a=sigma;
fs=115e3;
tic
parfor i=firstSignal:lastSignal
    fitted=fitMySignal(mistery_data,i,'no');
    sigma(i)=fitted.sigma;
    delta(i)=fitted.delta;
    t_c(i)=fitted.t_c;
    a(i)=fitted.a;
    % scale back
    a(i)=a(i)*max(abs(mistery_data{i})); % scale amplitude with max
    Ns=length(mistery_data{i});
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

% histogram 

figure()
%histogram(diam);
% to use M BIN
histogram(diam,50);
title('Amplitude distribution') % ?????
xlabel('Electric diameter [\mu m]')

%% SCATTER PLOT and LINEAR FIT
% scatter plot
scatter_fig=figure();
scatter(diam,shape)
xlabel('Electric diameter [\mu m]')
ylabel('Shape parameters')

% getline, su fig --> coordinate dei punti della poligonale 
% inpolygon, -->  
% %%
%%
L=40e-3; %[mm] L=40 [um]
velocity=L./delta; % where L [um] and delta [s]
figure();
scatter(diam,velocity)
 % ATTENZIONE MANC QUALCOSA NELLE UNIT° DI MISURA !!!!
xlabel('Electric diameter [\mu m]')
ylabel('Velocity [mm/s]')
%% 
[X,Y]=getline(scatter_fig);
[index_value]=inpolygon(diam,shape,X,Y);
selected_value=[diam(index_value),shape(index_value)];
fit_line=polyfit(selected_value(:,1),selected_value(:,2),1);
fitted_y=polyval(fit_line,selected_value(:,1));
hold on
plot(selected_value(:,1),fitted_y,'LineWidth',2)

%% 
% retta y=p1x+p2 -> x=y/p1 - p2/p1 -> - p2/p1 + y/p1
% c1+c2(shape)
p1=fit_line(1);
p2=fit_line(2);
c1=-(p2/p1);
c2=1/p1;

diam_corr=diam./(c1+c2*shape);
figure()
histogram(diam_corr,50);
figure()
scatter(diam_corr,shape)
xlim([0,2])

%% TEST

% retta y=p1x+p2 -> x=y/p1 - p2/p1 -> - p2/p1 + y/p1
% c1+c2(shape)
p1=fit_line(1);
p2=fit_line(2);
c1=-(p2/p1);
c2=1/p1;

diam_corr=diam./(c1+c2*shape);
figure()
histogram(diam_corr,50);
figure()
scatter(diam_corr,shape)
xlim([0,2])


%% 
n_bin_d=200;
n_bin_sig=200;
% limits from Errico_Caselli_SAB_2017 

electric_D_lim     =[4.5,10.5];
electric_D_norm_lim=[0.8,1.6];
velocity_lim       =[0.1,0.5];
sig_ov_del_lim     =[0.15,0.30];
figure()
histogram2(diam_corr,shape,'DisplayStyle','tile','ShowEmptyBins','on', ...
    'NumBins',[n_bin_d,n_bin_sig],'XBinLimits',electric_D_norm_lim,'YBinLimits',sig_ov_del_lim);
colorbar