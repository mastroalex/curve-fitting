% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Mastrofini Alessandro
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Medical Engineering - University of Rome Tor Vergata
% Physiological Systems Modeling and Simulation
% F. Caselli, MSSF A.Y. 2021/2022
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Curve fitting for impedance micro cytometers
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function fitted_signal=fitMySignal(mistery_data,index,WantPlot)
    arguments %check variable integrity
        mistery_data (:,:) 
        index (1,:) {mustBeNumeric,SizeCheck(index,mistery_data)} %check index is lower that array lenght (see later)
        WantPlot (1,:) char {mustBeMember(WantPlot,{'yes','no'})} = 'yes' % plot or not plot figure
    end

    fs=115e3; %[Hz] % sampling frequencies
    T_sampling=1/fs; 
       
    g = fittype('a*(exp(-((t-(t_c-delta/2)).^2/(2*sigma.^2)))-exp(-((t-(t_c+delta/2)).^2/(2*sigma.^2))))',...
                'independent','t',...
                'coefficient',{'sigma','delta','t_c','a'},'dependent','data_fit');
    
    % It's better to give initial guess
    a0=1;
    delta0=1/3;
    sigma0=1/10;
    t_c0=0.5;
    
    % need to normalize signal amplitude
    % need to rescale time to te maximum time (last element of the array)
    % this is to avoid numerical problem with small numbers
    data_fit=mistery_data{index};
    Ns=length(data_fit);
    data_fit_t=1e3*(0:Ns-1)/fs; 
    % normalize
    data_fit=data_fit/max(abs(data_fit));
    t=data_fit_t/data_fit_t(end);
    % setup options like method and initial guess
    opts=fitoptions('Method','NonLinearLeastSquare');
    opts.StartPoint=[sigma0 delta0 t_c0 a0];
    % fitting
    fitted_signal=fit(t',data_fit',g,opts); 
    if isequal(WantPlot,'yes') % plot fitting curve and sample
        figure()
        % scale back time and fitted_signal
        % plot(scaled time, scale signal)
        % to scale back fitted: first calc in t and than scale it
        plot(t.*data_fit_t(end),fitted_signal(t).*max(abs(data_fit)),'-r')
        hold on
        % plot scaled samples
        plot(t.*data_fit_t(end),data_fit.*max(abs(data_fit)),'*b')
        xlabel('time [ms]')
        title(['Fit for mistery\_data: ', 'signal n.',num2str(index)])
        legend({'fit','data'})
        ylabel('signal [{\mu}A]')
    end
end

function SizeCheck(index,mistery) % check that index value is lower than max number of signals
    % Test for equal size
    if index>max(size(mistery))
        eid = 'Size:notEqual';
        msg = 'Size of the index must be lower than mistery_data{} lenght.';
        throwAsCaller(MException(eid,msg))
    end
end
