clc
clear

%% 系统模块要求及其对应参数
% 1. 信道编解码     ('None', 'LDPC', 'Polar')
channel_coding_mode = 'None';
% 2. 交织	         ('None', 'Interleaving')
interleaving_mode = 'None';
% 3. 调幅调相	     该参数直接参与计算，只能写数字
QAM_mode = 64;
% 4. 调制解调系统，此项不更改	 ('OTFS', 'OFDM') 
Modem_mode_1 = 'OFDM';
Modem_mode_2 = 'OFDM';
% 5. 同步（时延）    默认存在，代码里暂时不考虑
synchronization_mode = 'default';
% 6. 信道估计 ('None','LS','LMMSE'，'embedded_pilot','ideal')
% 只做有导频的信道估计，LS，LMMSE为频域估计，即FT
% embedded_pilot是指OTFS特定有的DD域估计,即下面的导频必须在DD域
% 理想的信道估计算法中的导频插入依托LMMSE的块状导频实现
channel_estimation_mode_OTFS = 'embedded_pilot';
channel_estimation_mode_OFDM = 'LMMSE';
% 7. 导频所在域(时域估计DT，频域估计FT，时延多普勒域估计DD,频率多普勒FD域),
pilot_domain_OTFS = 'DD';
pilot_domain_OFDM = 'FT';
% 8. 导频格式(块状导频block，梳状导频comb_like，格状导频lattice，用于DD域OTFS估计的OTFS_delta导频)
pilot_form_OTFS= 'OTFS_delta';
pilot_form_OFDM= 'block';
% 9. 信道均衡 ('None','LS',迫零均衡算法'ZF', 'MMSE',以上信道均衡默认在FT域进行，即频域均衡，与导频所在域无关)
channel_equalization_mode_OTFS = 'LS';
channel_equalization_mode_OFDM = 'LS';
% 10. 保护间隔            ('None','CP','ZP','RCP','RZP')
interval_mode = 'None';
% 11.仿真所用信道模型 ('None','AWGN','Rayleigh','NTN_TDL','NTN_TDL_Moving','Multipath_delay','Singlpath_doppler','Multipath_doppler')
channel_module = 'NTN_TDL_Moving';


% 补充4：后续仿真为了方便统一，可能会选择OFDM、OTFS共存的模式，因此需要两个变量
% 补充6，7：LS，LMMSE均为频域估计，块状导频。未来可能还有别的模式，包括华为的均衡技术，时延多普勒域的均衡，等等

% 补充1.编码格式
if(strcmp(channel_coding_mode,'LDPC')) 
    LDPC_rate = 1/2;
    LDPC_codeword_length = 672;
elseif(strcmp(channel_coding_mode,'Polar')) 
    Polar_rate = 1/2;
    Polar_codeword_length = 256;
    Polar_info_length = Polar_codeword_length*Polar_rate;
    Polar_CRC_decode_length = 8;
end
% 补充6，9.信道估计域均衡
if strcmp(channel_estimation_mode_OTFS,'None')
    channel_equalization_mode_OTFS = 'None';
end
if strcmp(channel_estimation_mode_OFDM,'None')
    channel_equalization_mode_OFDM = 'None';
end
if strcmp(channel_estimation_mode_OTFS,'embedded_pilot')
    pilot_domain_OTFS = 'DD';
    pilot_form_OTFS= 'OTFS_delta';
end
if strcmp(channel_estimation_mode_OFDM,'embedded_pilot')
    pilot_domain_OFDM = 'DD';
    pilot_form_OFDM= 'OTFS_delta';
end
if strcmp(channel_estimation_mode_OTFS,'LS')||strcmp(channel_estimation_mode_OTFS,'LMMSE')
    pilot_domain_OTFS = 'FT';
end
if strcmp(channel_estimation_mode_OFDM,'LS')||strcmp(channel_estimation_mode_OFDM,'LMMSE')
    pilot_domain_OFDM = 'FT';
end
% 补充11.信道报错
if ~sum(strcmp(channel_module,{'None','AWGN','Rayleigh','NTN_TDL','NTN_TDL_Moving','Multipath_delay','Singlpath_doppler','Multipath_doppler'}))
    msg = '请选择信道环境（None, AWGN, Rayleigh, NTN_TDL, NTN_TDL_Moving, Multipath_delay, Singlpath_doppler, Multipath_doppler)';
    error(msg)
end

%% 数据帧格式参数 OTFS（OFDM同理）
% OTFS参数确定

% 1. N：单位帧的时隙数，列数，不包含前缀，包含导频
N = 64;
% 2. M：帧的载波数，行数
M = 64;
% 3. 子载波间隔：15k，单位Hz
delta_f = 15*10^3; 
% 4. 时隙间隔：1/15k，单位s
delta_t = 1/delta_f; 

% 时延抽头的尺度 % 单位 秒
one_delay_tap = 1/(M*delta_f);
% 多普勒抽头的尺度 % 单位 Hz
one_doppler_tap = 1/(N*delta_t);

%% 信道仿真参数
% 1.载波频率
car_frequency = 2*10^9; % S波段，单位Hz
% 2.卫星轨道高度
sat_altitude =600*10^3; % 单位 m
%卫星开始服务时的极限俯仰角
sat_edge_elevation_ang = 10;   %角度单位
sat_edge_elevation_rad = sat_edge_elevation_ang/360*(2*pi);%弧度单位


sat_ele_ang = 50;


% 地面UE最高移动速度
UE_max_speed = 500;     % 单位 km/h
% 信噪比
SNR_dB = 15:2.5:30;
SNR = 10.^(SNR_dB/10);
% 单次循环仿真的帧数
Simulation_frame = 1000;

%% 其它内容的初始化

% 一个波形承载的比特
QAM_bits = log2(QAM_mode);
% 一个波形符号的平均能量（这部分之前有错误），假定不同模式下，星座点离最近的点的距离是一致的，都为单位1
eng_sqrt = (QAM_mode==2)+(QAM_mode~=2)*sqrt((QAM_mode-1)/6*(2^2));
% 噪声的能量
sigma_2 = (abs(eng_sqrt)^2)./SNR;



% 计算消息码元以及导频序列在一帧中的位置(位置取决于前缀，后缀)%%%%%


% 前缀 后缀等的初始化
% 原文的加CP方式和5G不符，5G中CP不占用信息位，
% 在3Gpp OFDM系统中，存在normal CP和extemd CP，两种模式，extend CP只在60KHz子载波间隔才有，本代码不考虑
% 常规模式（normal CP）分为长短CP，短CP的长度占一个OFDM符号的7%，（长CP不考虑）
length_CP  = 0;
length_ZP  = 0;
length_RCP = 0;
length_RZP = 0;
if(strcmp(interval_mode,'ZP'))         
    length_ZP = ceil(M/16); % ZP length (required only for ZP-OTFS)
    length_interval = length_ZP;
elseif(strcmp(interval_mode,'CP'))
    length_CP = ceil(M/16); % CP length (required only for CP-OTFS) 
    length_interval = length_CP;
elseif strcmp(interval_mode,'RCP')
    length_RCP = ceil(M/16);
    length_interval = length_RCP;
elseif strcmp(interval_mode,'RZP')
    length_RZP = ceil(M/16);
    length_interval = length_RZP;
elseif strcmp(interval_mode,'None')
    length_interval = 0;
else
    msg = '选择OTFS的前缀模式: (None/ RZP / RCP / CP / ZP)';
    error(msg)
end

% 在时间——多普勒上加CP  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    length_CP2 = ceil(N/16);

% 每一帧能承载的符号（波形）数量，包含导频与数据，不包含前缀/后缀
symbols_perframe = M*N;
% 每一帧能承载的码元（比特）数量，包含导频与数据，不包含前缀
bits_perframe = symbols_perframe*QAM_bits;


% 信道估计与均衡用的导频初始化，以及除去导频后的剩余码元长度（导频会侵占数据位空间）
if strcmp(channel_estimation_mode_OTFS,'None')
    N_data_OTFS = N;
    M_data_OTFS = M;
else
    [M_pilot_OTFS,N_pilot_OTFS,M_data_OTFS,N_data_OTFS,pilot_Indx_OTFS,Data_Indx_OTFS] = function_pilot_initilization(M,N,channel_estimation_mode_OTFS,pilot_form_OTFS);
end
if strcmp(channel_estimation_mode_OFDM,'None')
    N_data_OFDM = N;
    M_data_OFDM = M;
else
    [M_pilot_OFDM,N_pilot_OFDM,M_data_OFDM,N_data_OFDM,pilot_Indx_OFDM,Data_Indx_OFDM] = function_pilot_initilization(M,N,channel_estimation_mode_OFDM,pilot_form_OFDM);
end

%%%% 这部分代码是有问题的，如果OFDM和OTFS占的导频空间不一致，请更改这部分的代码
if 1% (M_data_OTFS*N_data_OTFS ==M_data_OFDM*N_data_OFDM)
symbols_data_perframe = M_data_OTFS*N_data_OTFS;
bits_data_perframe = symbols_data_perframe*QAM_bits; % 留出导频空间
else
    msg = '如果OFDM和OTFS占的导频空间不一致，请更改导频初始化部分的代码';
    error(msg)
end


% FFT 矩阵的初始化（降低运算复杂度）
Fn_M_OTFS=dftmtx(M_data_OTFS);  % 生成DTF矩阵，矩阵大小为MxM，用于频率-时延的FFT转换，使用时左乘
Fn_M_OTFS=Fn_M_OTFS./norm(Fn_M_OTFS);  % 归一化的矩阵为酉矩阵，此时转秩等于逆
Fn_N_OTFS=dftmtx(N_data_OTFS);  % 生成DTF矩阵，矩阵大小为NxN，用于时间-多普勒的FFT转换，使用时右乘
Fn_N_OTFS=Fn_N_OTFS./norm(Fn_N_OTFS);  % 归一化的矩阵为酉矩阵，此时转秩等于逆

Fn_M_OFDM=dftmtx(M_data_OFDM);  % 生成DTF矩阵，矩阵大小为MxM，用于频率-时延的FFT转换，使用时左乘
Fn_M_OFDM=Fn_M_OFDM./norm(Fn_M_OFDM);  % 归一化的矩阵为酉矩阵，此时转秩等于逆
Fn_N_OFDM=dftmtx(N_data_OFDM);  % 生成DTF矩阵，矩阵大小为NxN，用于时间-多普勒的FFT转换，使用时右乘
Fn_N_OFDM=Fn_N_OFDM./norm(Fn_N_OFDM);  % 归一化的矩阵为酉矩阵，此时转秩等于逆



Fn_M_pilot = dftmtx(M); % 用于有导频信息的频率-时延域转换，使用时左乘
Fn_M_pilot=Fn_M_pilot./norm(Fn_M_pilot);
Fn_N_pilot = dftmtx(N); % 用于有导频信息的时间-多普勒域转换，使用时右乘
Fn_N_pilot=Fn_N_pilot./norm(Fn_N_pilot);



% 计算一帧之中数据码元数量，取决于导频和编解码

if strcmp(channel_coding_mode,'None')
    info_bits_length = bits_data_perframe;
elseif strcmp(channel_coding_mode,'LDPC')
    % LDPC_rate =  ;                % LDPC数据码元占比，上面有，原文只提供1/2LDPC编码，3/4的LDPC编码没见到
    % LDPC_codeword_length = ;      % LDPC编码后的一个码块的码字长度，上面有，可选(672 / 3840)，其他的需要自己编写H矩阵
    LDPC_info_length  = LDPC_codeword_length*LDPC_rate;              % 一个LDPC码块的数据码元数量
    LDPC_trans_blocks = floor(bits_data_perframe/LDPC_codeword_length);  % 一个数据帧能够承载的LDPC码块数，只传整数，小数部分不传。
    LDPC_bits_length  = LDPC_trans_blocks*LDPC_codeword_length;      % 多个LDPC码块填满一帧，码块所占的比特数（因为码块不能整除数据帧，所以二者并不相等）
    info_bits_length  = LDPC_trans_blocks*LDPC_info_length;          % 多个LDPC码块填满一帧，数据码元所占一帧的比特数，即码块长度乘数据码元所占百分比
elseif strcmp(channel_coding_mode,'Polar')    
    Polar_info_length  = Polar_codeword_length*Polar_rate;              % 一个Polar码块的数据码元数量
    Polar_trans_blocks = floor(bits_data_perframe/Polar_codeword_length);  % 一个数据帧能够承载的Polar码块数，只传整数，小数部分不传。
    Polar_bits_length  = Polar_trans_blocks*Polar_codeword_length;      % 多个Polar码块填满一帧，码块所占的比特数（因为码块不能整除数据帧，所以二者并不相等）
    info_bits_length   = Polar_trans_blocks*Polar_info_length;          % 多个Polar码块填满一帧，数据码元所占一帧的比特数，即码块长度乘数据码元所占百分比
else
    msg = '选择编码方式 （None, LDPC, Polar）';
    error(msg)
end

%对编解码矩阵的初始化
if strcmp(channel_coding_mode,'None')
    % 无需初始化
elseif strcmp(channel_coding_mode,'LDPC')%%%%%%%%%%%%%%%%%%%%%%%%% H矩阵生成，目前暂时用的原代码，这部分要改
    [hEnc,hDec,hDec_coded_soft,hDec_coded_hard]=LDPC_system_objects(LDPC_rate,LDPC_codeword_length); 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(channel_coding_mode,'Polar')    
    % 无
end


% 误码计数的初始化
if strcmp(Modem_mode_1,'OTFS')||strcmp(Modem_mode_2,'OTFS')
    result_err_ber_OTFS = zeros(1,length(SNR_dB));% 每一帧的误码计数
    result_avg_ber_OTFS = zeros(1,length(SNR_dB));% 每一次循环的误码平均值
    result_err_fer_OTFS = zeros(1,length(SNR_dB));% 误帧
    result_avg_fer_OTFS = zeros(1,length(SNR_dB));
end

if strcmp(Modem_mode_1,'OFDM')||strcmp(Modem_mode_2,'OFDM')
    result_err_ber_OFDM = zeros(1,length(SNR_dB));% 每一帧的误码计数
    result_avg_ber_OFDM = zeros(1,length(SNR_dB));% 每一次循环的误码平均值
    result_err_fer_OFDM = zeros(1,length(SNR_dB));% 误帧
    result_avg_fer_OFDM = zeros(1,length(SNR_dB));
end

%代码运行时间的计数初始化


% 仿真次数计数初始化 
simulation_frame = zeros(1,length(SNR_dB));

% 嵌入式导频的相对于信号的幅度
embeded_pilot_amp = sqrt(N);

% N为64的话，取sqrt（N），则导频能量和符号能量差了N倍，即18dB
% N为64的话，取N，则导频能量和符号能量差了N*N倍，即36dB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 仿真参数初始化完成 %%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 仿真开始 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 信号生成，列向量，长度不固定
    info_tx_bits = randi([0,1],info_bits_length,1); % 数据长度取决于所用编码，以及导频
    % info_tx_bits = ones(info_bits_length,1); % 数据长度取决于所用编码，以及导频
    % info_tx_bits = zeros(info_bits_length,1); % 数据长度取决于所用编码，以及导频
% 编码：无编码，输出尺寸（M*N_data*QAM，1）
        info_tx_bits_coded = info_tx_bits;    % 尺寸为M_data*N_data*QAM
% 交织：无交织，输出尺寸（M_data*N_data*QAM，1）
        info_tx_bits_interleaved = info_tx_bits_coded;
% 调幅调相（QAM）输出尺寸（1, M_data*N_data）
        info_tx_symbols_QAM=qammod(reshape(info_tx_bits_interleaved,QAM_bits,symbols_data_perframe), QAM_mode,'gray','InputType','bit');  %格雷码  
% 串并转换 sreializer to parallel ,输出尺寸 (M_data,N_data)
        info_tx_symbols_S2P_OTFS = reshape(info_tx_symbols_QAM,M_data_OTFS,N_data_OTFS); % 按列读取，按列摆放
% 数字调制（OTFS/OFDM）,输出尺寸 (M_data,N_data)
        info_tx_OTFS_DT = info_tx_symbols_S2P_OTFS*Fn_N_OTFS';% DT代表delay_time域 
      % info_tx_OTFS_DT = Fn_M_OTFS'*info_tx_symbols_S2P_OTFS;
% 插入导频（频域导频），输出尺寸 (M,N)    
        pilot_symbols_OTFS = function_pilot_generation(QAM_mode,M_data_OTFS,N_data_OTFS,pilot_form_OTFS,embeded_pilot_amp);
        if strcmp(channel_estimation_mode_OTFS,'None')
            info_tx_OTFS_pilot_DT = info_tx_OTFS_DT;
        else
            info_tx_OTFS_pilot_DT = function_pilot_insert(M,N,Fn_M_OTFS,Fn_N_OTFS,Fn_M_pilot,Fn_N_pilot,info_tx_OTFS_DT,pilot_symbols_OTFS,pilot_domain_OTFS,pilot_form_OTFS,Data_Indx_OTFS,pilot_Indx_OTFS);
        end
% 插入保护间隔（循环前缀），CP加在前面，ZP加在后面
        info_tx_OTFS_interval_DT = info_tx_OTFS_pilot_DT;
        info_tx_OTFS_RCP_DT  = [];
        info_tx_OTFS_RZP_DT  = [];
% 并串转换，parallel to sreializer ,输出尺寸 (MxN,1),或者((M+CP/ZP)xN,1),或者(MxN+RCP/RZP,1)
        info_tx_OTFS_P2S_DT = [info_tx_OTFS_RCP_DT; reshape(info_tx_OTFS_interval_DT,numel(info_tx_OTFS_interval_DT),1); info_tx_OTFS_RZP_DT]; 




%% 仿真中段（信道部分）
% 信道环境生成


        [T_rms,serv_t_percent,channel_amplitude_coef_taps,channel_delay_taps,channel_doppler_taps,channel_taps]=...
            function_Gen_NTN_TDL_series_channel_parameters_DD(M,N,car_frequency,delta_f,delta_t,UE_max_speed,sat_ele_ang, channel_module);
%         channel_delay_taps    = [0 1 2];
%         channel_doppler_taps  = [1 1 1];

% 信号同步 抵消掉大时延造成的影响，具体实现，大延迟靠计算提前量（TA估计），小延迟由导频完成，这里进行了简化。如有需要另行编写
    if strcmp(synchronization_mode, 'default')
        channel_delay_taps = channel_delay_taps - channel_delay_taps(1);
    end
    % channel_delay_taps = [0,1,2,3];
% 信号同步 抵消掉大多普勒造成的影响，这里进行了简化。如有需要另行编写
%     if strcmp(synchronization_mode, 'default')
%         channel_doppler_taps = (channel_doppler_taps - channel_doppler_taps(1));
%     end
% 合并相同时延的参数（不在实验项目书的步骤中）
    channel_delay_set=unique(channel_delay_taps); 
    channel_delay_max=max(channel_delay_set);



% 生成信道的时域表示i_SNR
    Channel_T=function_Gen_delay_time_channel...
    (M,N,channel_delay_taps,channel_doppler_taps,channel_amplitude_coef_taps,length_interval,interval_mode); % delay time域

%     Channle_F = fft(Channel_T)/N;
%     a = find(Channle_F>1);
%     Channel_DT = reshape(Channel_T,size(info_tx_OTFS_interval_DT));
%     Channel_DD = Channel_DT*Fn_N_pilot;
%     Channel_FD = Fn_M_pilot*Channel_DD;

    Channel_HDT=function_Gen_delay_time_channelH...
    (M,N,channel_delay_taps,channel_doppler_taps,channel_amplitude_coef_taps,length_interval,interval_mode); % 用于理想信道估计的信道矩阵

    Channel_HDD = Channel_HDT*Fn_N_pilot;

% 信号通过信道
    i_SNR = 1;
    if strcmp(channel_module,'None')
            info_channel_out_OTFS_DT = info_tx_OTFS_P2S_DT; 
    else
            info_channel_out_OTFS_DT = ...
                function_Gen_channel_output_DT(channel_delay_set,Channel_T,info_tx_OTFS_P2S_DT,sigma_2(i_SNR)); 
    end



%% 仿真末端（接收部分）

% 高频解调（无）
% 低频解调（无）
% 信号同步（小时延） 上文在大时延处 已经进行了相应的简化，此处也不用编写 


% 去掉保护间隔/循环前缀，同时串并转换，这里默认知道矩阵尺寸, 输出尺寸（M，N）

            info_rx_OTFS_S2P_DT  = reshape(info_channel_out_OTFS_DT,size(info_tx_OTFS_interval_DT));
            info_rx_OTFS_deinterval_DT = info_rx_OTFS_S2P_DT;

        info_rx_FT= Fn_M_pilot*info_rx_OTFS_deinterval_DT;




% 信道估计
    % channel_estimation_mode_OTFS ='None';
    if strcmp(channel_estimation_mode_OTFS,'None')
        info_rx_OTFS_estimation_DT = info_rx_OTFS_deinterval_DT; 
    else
        [info_rx_OTFS_domain,channel_H_est_OTFS_domain] = ...
            function_channel_estimation...
            (M_pilot_OTFS,N_pilot_OTFS,Fn_M_pilot,Fn_N_pilot,...
            info_rx_OTFS_deinterval_DT,pilot_symbols_OTFS,...
            QAM_mode,pilot_form_OTFS,pilot_domain_OTFS,channel_estimation_mode_OTFS,...
            pilot_Indx_OTFS,SNR(i_SNR),sigma_2(i_SNR), Channel_HDT,interval_mode,length_interval);
    end





% 信道均衡,输出尺寸（M_data，N_data）
    if strcmp(channel_estimation_mode_OTFS,'None')
        info_rx_OTFS_equalization_DT = info_rx_OTFS_estimation_DT;
    else
        info_rx_OTFS_equalization_DT = function_channel_equalization...
            (Fn_M_OTFS,Fn_N_OTFS,Fn_M_pilot,Fn_N_pilot,...
            channel_equalization_mode_OTFS,pilot_domain_OTFS,pilot_form_OTFS,...
            info_rx_OTFS_domain,channel_H_est_OTFS_domain,...
            Data_Indx_OTFS,sigma_2(i_SNR));
    end

% 数字解调（OTFS/OFDM）,输出尺寸（M，N_data）

        info_rx_deOTFS_DD = info_rx_OTFS_equalization_DT * Fn_N_OTFS; % 归一化的FFT矩阵是酉矩阵，共轭转秩等于矩阵的逆


% 并串转换    输出尺寸（1,M*N_data）
        info_rx_OTFS_P2S = reshape(info_rx_deOTFS_DD,1,numel(info_rx_deOTFS_DD));



% 解QAM，输出尺寸（QAM*M*N_data,1）
if strcmp(channel_coding_mode,'None')
    if exist('info_tx_OTFS_DT','var') 
        info_rx_bits_OTFS_deQAM = qamdemod(info_rx_OTFS_P2S,QAM_mode,'gray','OutputType','bit');
        info_rx_bits_OTFS_deQAM = reshape(info_rx_bits_OTFS_deQAM,numel(info_rx_bits_OTFS_deQAM),1);
    end
else
    if exist('info_tx_OTFS_DT','var') 
        info_rx_bits_OTFS_deQAM = qamdemod(info_rx_OTFS_P2S,QAM_mode,'gray','OutputType','llr');
        info_rx_bits_OTFS_deQAM = reshape(info_rx_bits_OTFS_deQAM,numel(info_rx_bits_OTFS_deQAM),1);
    end
end


% 解交织，输出尺寸（QAM*M*N_data,1）
    if strcmp(interleaving_mode, 'Interleaving')
            info_rx_bits_OTFS_deinterleaved = info_rx_bits_OTFS_deQAM(RND_Intrlv_rev); 
    else
            info_rx_bits_OTFS_deinterleaved = info_rx_bits_OTFS_deQAM; 
    end


% 解码
            info_rx_bits_OTFS_decoded = info_rx_bits_OTFS_deinterleaved;    % 尺寸为MxNxQAM





scatterplot(info_tx_symbols_QAM); 
hold on;
axis([-8 8 -8 8])
grid on
title('发送端信号星座图','Interpreter','none')
hold off

scatterplot(info_rx_OTFS_P2S); 
hold on;
axis([-8 8 -8 8])
grid on
title([channel_module,'环境下OTFS系统接收端信号星座图，信噪比',num2str(SNR_dB(i_SNR)),'dB'],'Interpreter','none')
