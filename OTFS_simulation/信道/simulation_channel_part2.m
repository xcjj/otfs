

% 生成信道的时域表示
    Channel_DT=function_Gen_delay_time_channel...
    (M,N,channel_delay_taps,channel_doppler_taps,channel_amplitude_coef_taps,length_interval,interval_mode); % delay time域
    
    Channel_HDT=function_Gen_delay_time_channelH...
    (M,N,channel_delay_taps,channel_doppler_taps,channel_amplitude_coef_taps,length_interval,interval_mode); % 用于理想信道估计的信道矩阵

% 信号通过信道
if strcmp(channel_module,'None')
    if exist('info_tx_OTFS_DT','var') 
        info_channel_out_OTFS_DT = info_tx_OTFS_P2S_DT; 
    end
    if exist('info_tx_OFDM_DT','var')
        info_channel_out_OFDM_DT = info_tx_OFDM_P2S_DT; 
    end
else
    if exist('info_tx_OTFS_DT','var') 
        info_channel_out_OTFS_DT = function_Gen_channel_output_DT(channel_delay_set,Channel_DT,info_tx_OTFS_P2S_DT,sigma_2(i_SNR)); 
    end
    if exist('info_tx_OFDM_DT','var')
        info_channel_out_OFDM_DT = function_Gen_channel_output_DT(channel_delay_set,Channel_DT,info_tx_OFDM_P2S_DT,sigma_2(i_SNR)); 
    end
end
