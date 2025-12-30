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
SNR_dB = 5:2.5:40;
SNR = 10.^(SNR_dB/10);
% 单次循环仿真的帧数
% Simulation_frame = 1000;
simulation_frame = zeros(1,length(SNR_dB));
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
%     length_CP2 = ceil(N/16);

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
    % % 导频格式(块状导频block，格状导频lattice，梳状导频comb_like)
%     if strcmp(pilot_form_OTFS,'block')
%         pilot_Inter=8; % 导频间隔 pilot interval,导频间隔是8意味着导频间计数需要+9
%         % 计算导频和数据数量
%         M_pilot = 0;
%         N_pilot=ceil(N/(pilot_Inter+1)); % 向上取整
%         if rem(N,pilot_Inter+1)==0 % 求余运算，如果被整除，则在末尾加一个导频
%             N_pilot=N_pilot+1;
%         end
%         % 数据码元的数量
%         M_data = M;
%         N_data = N-N_pilot;
%         % 导频与信息在一帧中位置的计算
%         pilot_Indx=zeros(1,N_pilot);
%         for i=1:N_pilot-1
%             pilot_Indx(1,i)=(i-1)*(pilot_Inter+1)+1;
%         end
%         pilot_Indx(1,N_pilot)=N; % 默认最后一位必须是导频
% 
%         Data_Indx = 1:N;
%         Data_Indx(pilot_Indx) =[];
%     elseif strcmp(pilot_form_OTFS,'lattice')
%         pilot_Inter=8; % 导频间隔 pilot interval,导频间隔是8意味着导频间计数需要+9
%         % 计算导频和数据数量
%         M_pilot=ceil(M/(pilot_Inter+1)); % 向上取整
%         if rem(M,pilot_Inter+1)==0 % 求余运算，如果被整除，则在末尾加一个导频
%             M_pilot=M_pilot+1;
%         end
%         N_pilot = 0;
%         % 数据码元的数量
%         N_data = N;
%         M_data = M-M_pilot;
%         % 导频与信息在一帧中位置的计算
%         pilot_Indx=zeros(M_pilot,1);
%         for i=1:M_pilot-1
%             pilot_Indx(i,1)=(i-1)*(pilot_Inter+1)+1;
%         end
%         pilot_Indx(M_pilot,1)=M; % 默认最后一位必须是导频
%         Data_Indx = (1:M)';
%         Data_Indx(pilot_Indx) =[];
%     elseif strcmp(pilot_form_OTFS,'comb_like')
%         %%%%%%
%     end

end
if strcmp(channel_estimation_mode_OFDM,'None')
    N_data_OFDM = N;
    M_data_OFDM = M;
else
    [M_pilot_OFDM,N_pilot_OFDM,M_data_OFDM,N_data_OFDM,pilot_Indx_OFDM,Data_Indx_OFDM] = function_pilot_initilization(M,N,channel_estimation_mode_OFDM,pilot_form_OFDM);
end

%%%% 这部分代码是有问题的，如果OFDM和OTFS占的导频空间不一致，请更改这部分的代码
if (M_data_OTFS*N_data_OTFS ==M_data_OFDM*N_data_OFDM)
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


% 嵌入式导频的相对于信号的幅度
embeded_pilot_amp = N;
% N为64的话，取sqrt（N），则导频能量和符号能量差了N倍，即18dB
% N为64的话，取N，则导频能量和符号能量差了N*N倍，即36dB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 仿真参数初始化完成 %%%%%%%%%%%%%%%%%%%%%%%%%%%

