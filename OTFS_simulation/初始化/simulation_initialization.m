%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 仿真参数初始化 %%%%%%%%%%%%%%%%%%%%%%%%%%%

% 截止到20230524，本代码更新如下：

% 加入polar码编码，使用matlab自带编码，注意：Polar码解码用时会很长

% 信道估计可以在NTN——TDL环境，以及地面高速移动时的NTN——TDL环境下使用
% LS,LMMSE为频域估计，OTFS_delta为时延/多普勒估计，选项7导频所在域会根据选项6自动更改
% OTFS_delta信道估计的相关文献附在压缩包中，该估计对噪声非常敏感，低信噪比情况下几乎不可用

% OTFS与OFDM的信道估计与均衡可以独立设置参数
% LS，LMMSE的信道估计的导频格式分为块状，格状，梳状。其中梳状导频没写，暂时用不了，OTFS的信道估计有自己的导频格式

% 信道均衡都是频域均衡，分为ZF，LS，MMSE，区别如下，正常选择LS均衡即可，（MMSE均衡也没做好）
% https://blog.csdn.net/DqiangLiu/article/details/127416035

% 加入理想的信道估计结果

% 由于更改较大，部分代码跑起来可能存在错误




%% 系统模块要求及其对应参数
% 1. 信道编解码     ('None', 'LDPC', 'Polar')
channel_coding_mode = 'LDPC';
% 2. 交织	         ('None', 'Interleaving')
interleaving_mode = 'Interleaving';
% 3. 调幅调相	     该参数直接参与计算，只能写数字
QAM_mode = 64;
% 4. 调制解调系统，此项不更改	 ('OTFS', 'OFDM') 
Modem_mode_1 = 'OTFS';
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
interval_mode = 'CP';
% 11.仿真所用信道模型 ('None','AWGN','Rayleigh','NTN_TDL','NTN_TDL_Moving','Multipath_delay','Singlpath_doppler','Multipath_doppler')
channel_module = 'NTN_TDL_Moving';


% 补充4：后续仿真为了方便统一，可能会选择二者共存的模式，因此需要两个变量
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