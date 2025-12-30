%% 仿真中段（信道部分）
% 信道环境生成
    if strcmp(channel_module,'AWGN')||strcmp(channel_module,'None')
        channel_delay_taps = 0; 
        channel_taps = length(channel_delay_taps); 
        channel_amplitude_coef_taps = 1;
        channel_doppler_taps = 0;
    elseif strcmp(channel_module,'Rayleigh')    % 函数：rayleighchan 可以产生时变的瑞利信道，这里没用
        channel_delay_taps = 0; 
        channel_taps = length(channel_delay_taps); 
        channel_amplitude_coef_taps = 1*(sqrt(1/2) * (randn(1,channel_taps)+1i*randn(1,channel_taps)));
        channel_doppler_taps = 0;
    elseif sum(strcmp(channel_module,{'NTN_TDL','NTN_TDL_Moving','Multipath_delay','Singlpath_doppler','Multipath_doppler'}))
        [T_rms,serv_t_percent,channel_amplitude_coef_taps,channel_delay_taps,channel_doppler_taps,channel_taps]=function_Gen_NTN_TDL_series_channel_parameters_DD(M,N,car_frequency,delta_f,delta_t,UE_max_speed,sat_ele_ang, channel_module);
        % channel_delay_taps = [0 0 0 1];
    end

% 信号同步 抵消掉大时延造成的影响，具体实现，大延迟靠计算提前量（TA估计），小延迟由导频完成，这里进行了简化。如有需要另行编写
    if strcmp(synchronization_mode, 'default')
        channel_delay_taps = channel_delay_taps - channel_delay_taps(1);
    end
% 信号同步 抵消掉大多普勒造成的影响，这里进行了简化。如有需要另行编写
    if strcmp(synchronization_mode, 'default')
        channel_doppler_taps = channel_doppler_taps - channel_doppler_taps(1);
    end

% 合并相同时延的参数（不在实验项目书的步骤中）
    channel_delay_set=unique(channel_delay_taps); 
    channel_delay_max=max(channel_delay_set);

