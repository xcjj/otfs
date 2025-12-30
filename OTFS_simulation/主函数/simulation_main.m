clc
clear
rng('shuffle')
% rng是指定随机数生成器的规则或种子，其中rng('shuffle')是指根据当前时间初始化随机数生成器。

% 仿真系统的更改在 simulation_initialization.m 文件中
% 增加仿真进程显示
% 现在代码循环次数也进行了更改，低信噪比环境下只会执行200帧，每一轮仿真帧数随着信噪比增高而增高


% a_system_simulation_main.m为整体的文件，也可以用它跑仿真，内容是一样的
% fer_count = zeros(16,length(SNR_dB));

i=0;
simulation_initialization;
% embedded_pilot是指OTFS特定有的DD域估计,即下面的导频必须在DD域
channel_estimation_mode_OTFS = 'embedded_pilot';
channel_estimation_mode_OFDM = 'LMMSE';
% 7. 导频所在域(时域估计DT，频域估计FT，时延多普勒域估计DD,频率多普勒FD域),
pilot_domain_OTFS = 'DD';
pilot_domain_OFDM = 'FT';
% 8. 导频格式(块状导频block，梳状导频comb_like，格状导频lattice，用于DD域OTFS估计的OTFS_delta导频)
pilot_form_OTFS= 'OTFS_delta';
pilot_form_OFDM= 'block';
simulation_initialization_others;



for i_SNR = 1:length(SNR_dB)
    for i_frame = 1:10000 %: Simulation_frame
        simulation_tx;
        simulation_channel_part1;
        simulation_channel_part2;
        simulation_rx;


        simulation_error_count;


        % 进程显示
        process_display;
        
        if ~xor(exist('info_tx_OTFS_DT','var'),(result_err_fer_OTFS(1,i_SNR)>=200))&& ~xor(exist('info_tx_OFDM_DT','var'),(result_err_fer_OFDM(1,i_SNR)>=200))
            break;
        end

    end
    
end

i = i+1;
ber_count_OTFS(i,:) = result_avg_ber_OTFS; 
fer_count_OTFS(i,:) = result_avg_fer_OTFS;
ber_count_OFDM(i,:) = result_avg_ber_OFDM; 
fer_count_OFDM(i,:) = result_avg_fer_OFDM;

%  simulation_figure;

