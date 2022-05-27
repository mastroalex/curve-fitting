%% Lezione 6 - 16/03/2022

load mistery.mat

%% Step 1 visulizzare i segnali
close all
figure()
hold on
starting_signal=1;
ending_signal=5;
for i=starting_signal:ending_signal
    plot(mistery_data{i},'*')
end
%% Step 2 plottare il template
t_c=0; % given by prof.ssa 
a=2; % see rference Article
sigma=2;
delta=8;

% define bipolar gaussian
f= @(t) a*(exp(-((t-(t_c-delta/2)).^2/(2*sigma.^2)))-exp(-((t-(t_c+delta/2)).^2/(2*sigma.^2))));

figure()
t=-100:0.05:100;
plot(t,f(t))
xlim([-20;20])
title('Bipolar gaussian example')
%% Considerare l'intervallo temporale [-10;10]

figure()
t=-10:0.01:10;
plot(t,f(t))
title('Bipolar gaussian example')

%%
% freqeunza di campionamento 115 KHz
fs=115e3; %[Hz]
T_sampling=1/fs; 

index=1;
data_sampled=mistery_data{index}; %consider only the signal indicized by index
Ns=length(data_sampled); 
time_data=1e3*(0:Ns-1)/fs; % time data. In in ms (*1e3)
figure()
plot(time_data,data_sampled)
xlabel('[ms]')
title(['mistery\_data: ', 'signal n.',num2str(index)])
% distanza in sensing zone a distanza di 40 micron 
% velocità 40 micron/ tempo picco-picco
%% STUDIARE I COMANDI fittype e fit

g = fittype('a*(exp(-((t-(t_c-delta/2)).^2/(2*sigma.^2)))-exp(-((t-(t_c+delta/2)).^2/(2*sigma.^2))))',...
            'independent','t',...
            'coefficient',{'sigma','delta','t_c','a'},'dependent','data_fit');

% It's better to give initial guess
% suggested by prof
a0=1;
delta0=1/3;
sigma0=1/10;
t_c0=0.5;

% need to normalize signal amplitude
% need to rescale time to te maximum time (last element of the array)
index=1;
data_fit=mistery_data{index};
Ns=length(data_fit);
data_fit_t=1e3*(0:Ns-1)/fs; 
% normalizzo
data_fit=data_fit/max(abs(data_fit));
t=data_fit_t/data_fit_t(end);
opts=fitoptions('Method','NonLinearLeastSquare');
opts.StartPoint=[sigma0 delta0 t_c0 a0];

%%
fitted=fit(t',data_fit',g,opts); % aggiungere fitoptions NonLineaLeastSquar 
% dentro fitted troviamo anche i parametri con le confidenze

%%
figure()
plot(fitted)
hold on
plot(t,data_fit,'*b')
xlabel('time [ms]')
title(['Fit for mistery\_data: ', 'signal n.',num2str(index)])
legend({'fit','data'})
ylabel('signal [A]')

%% Created a function to fit index signal and plot it vs fitted function 
% use like fitMySignal(index,mistery_data,WantPlot)
% WantPlot to 'no' and no plot
% WantPlot to 'yes' return also plot figure
% default is set to 'yes' to obtaine plotted figure
% the function return plot figure and fitting output
index=1000000;

%fitted=fitMySignal(mistery_data,index,'no');
fitted=fitMySignal(mistery_data,index);