% 原有模型的时延频偏只存在整数部分，现针对分数部分频偏进行建模
% 卫星通信信道仿真，卫星参数，场景参数：3Gpp波束凝视模型
% 参数需求：
%     轨道高度：600km
%     卫星开始服务的俯仰角度：10度
%     接收机底噪 7db(暂时无用)
% 自行设定的参数（具体值由函数外部确定）
%     Δf = 15k
%     Δt = 1/Δf
%     载波数M = 64
%     时隙数N = 64
%     调幅调相= 64QAM
%     所用频段 2GHz     
%     接收机信噪比：（取决于外部变量）
% 可计算参数：
%     卫星速度
%     过顶时间
%     多普勒频移（实时）
% 
% 
% 对时变信道的建模
% 
% NTN—TDL—D模型，（38.811 6.9.2，P91）
% 目前阶段参数缺失，仅支持50度俯仰角的相关计算。
% 2022.11.22 董开原
% 

% 原代码来源
%  [R1]. T. Thaj and E. Viterbo, "Low Complexity Iterative Rake Decision Feedback Equalizer for Zero-Padded OTFS Systems," in IEEE Transactions on Vehicular Technology, vol. 69, no. 12, pp. 15606-15622, Dec. 2020, doi: 10.1109/TVT.2020.3044276.
%  [R2]. T. Thaj and E. Viterbo,``Low Complexity Iterative Rake Detector for Orthogonal Time Frequency Space Modulation拻 2020 IEEE Wireless Communications and Networking Conference (WCNC), 2020, pp. 1-6, doi: 10.1109/WCNC45663.2020.9120526.
%  [R3]. Y. Hong, T. Thaj, E. Viterbo, ``Delay-Doppler Communications: Principles and Applications'', Academic Press, 2022, ISBN:9780323850285





function [T_rms,serv_t_percent,chan_amplitude_coef_taps,delay_taps,doppler_taps,taps]=function_Gen_NTN_TDL_series_channel_parameters_DD(M,N,car_frequency,delta_f,delta_t,UE_max_speed,sat_ele_ang, channel_module)

M_earth = 5.9742*10^24;       % 地球质量kg
R_earth = 6371*10^3;          % 地球半径m
G = 6.67259*10^-11;           % 万有引力常量N·m^2/kg^2
c = 299792458;                % 光速 m/s


% car_frequency = 2*10^9; % 载波频率
% delta_f = 15*10^3;      % 子载波间隔 Hz
% delta_t = 1/delta_f;    % 时隙长度
% N = 64;                 % 单位帧的时隙数
% M = 64;                 % 帧的载波数
% max_speed = 500;        % 地面移动速度单位km/h
% channel_module = 'NTN_TDL';
% sat_ele_ang = 50;

sat_alt =600*10^3;      % 卫星轨道高度 satellite orbital altitude(m)

%卫星服务时的极限俯仰角 Satellite elevation angle
sat_edg_ele_ang = 10;   %角度单位
sat_edg_ele_rad = sat_edg_ele_ang/360*(2*pi);%弧度单位




%%%%%%%%%%%%%%%%%%%%%%%%%%相关参数计算部分开始%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%卫星速度 satellite velocity m/s
sat_v= sqrt(G*M_earth/(R_earth+sat_alt));

%计算最大俯仰角的情况下卫星-地心-用户的角度(地心角 geocentric angle)
syms sina cosa;
f1 = ((sat_alt+R_earth)*cosa-R_earth)-tan(sat_edg_ele_rad)*(sat_alt+R_earth)*sina;
f2 = 1-cosa^2-sina^2;
[sina,cosa] = solve(f1,f2,sina,cosa);
geocent_max_rad = double(asin(sina(2)));%弧度单位
%geocent_max_angle = geocent_max_rad/(2*pi)*360;%角度单位

%卫星绕地一周所需时间 Time required for satellite to circle the earth 
sat_cir_time = 2*pi*(R_earth+sat_alt)/sat_v;%单位s
% 卫星服务UE最大时间；Maximum service time
Max_serv_t = sat_cir_time *2*geocent_max_rad/(2*pi);
%地心角变化速度, Geocentric angular velocity（匀速）
geocent_v_rad = 2*pi/sat_cir_time; % 单位 rad/s
%geocent_v_ang = geocent_v_rad/(2*pi)*360;



%%%%%%%%%%%%%%%%%%%至此相关参数计算部分完成,信道参数生成%%%%%%%%%%%%%%%%%%%

% 已知服务时间求相关参数
%服务时间，从-10度俯仰角接入开始，到+10度俯仰角接入结束，service time，取值范围0~Max_serv_t
% %serv_t = 0;%该参数可以手动输入
% serv_t = floor(serv_t_percent/100*Max_serv_t);%该参数为外部输入
% sat_geocent_rad = -geocent_max_rad+serv_t*geocent_ang_v;   %已知服务时间求地心角，弧度单位
% sat_geocent_ang = sat_geocent_rad/(2*pi)*360;
% %卫星与用户之间的距离 Distance between Satellite多径 多径and UE 
% dis_sat_UE= (sat_alt+R)*sin(sat_geocent_rad)/cos(sat_ele_rad);%单位m
% %此时的俯仰角计算
% sat_ele_rad =  acos(-(sat_alt+R_earth)*sin(sat_geocent_rad)/dis_sat_UE);
% sat_ele_ang = sat_ele_rad/(2*pi)*360;
% % % 几种距离的计算方式
% % dis_sat_UE= (sat_alt+R_earth)*sin(sat_geocent_rad)/cos(sat_ele_rad);%单位m
% % dis_sat_UE= sqrt(R_earth^2*sin(sat_ele_rad)^2+sat_alt^2+2*sat_alt*R_earth)-R_earth*sin(sat_ele_rad);
% % dis_sat_UE = sqrt(((sat_alt+R_earth)*sin(sat_geocent_rad))^2+((sat_alt+R_earth)*cos(sat_geocent_rad)-R_earth)^2);%单位m

% 已知俯仰角求相关参数 sat_ele_ang(satellite elevation angle)
sat_ele_rad = sat_ele_ang/360*(2*pi);% 俯仰角的弧度单位
% 卫星与用户之间的距离 Distance between Satellite and UE
dis_sat_UE= sqrt(R_earth^2*sin(sat_ele_rad)^2+sat_alt^2+2*sat_alt*R_earth)-R_earth*sin(sat_ele_rad);% 单位m
% 卫星与UE的地心角
sat_geocent_rad = asin(-cos(sat_ele_rad)*dis_sat_UE/(sat_alt+R_earth)) ; % 弧度单位
%sat_geocent_ang = sat_geocent_rad/(2*pi)*360;% 角度单位
% 此时的服务时间占比（正过顶情况下）
serv_t = (sat_geocent_rad + geocent_max_rad)/geocent_v_rad;   %已知服务时间求地心角，弧度单位
serv_t_percent = floor(serv_t/Max_serv_t*10000)/100; %精度小数点后两位
%多普勒频移计算
doppler_sat = -car_frequency/c*sat_v*cos(sat_ele_rad+sat_geocent_rad); %帧头处的多普勒偏移
%时延参数的计算(传播时延)
delay_sat = dis_sat_UE/c; %单位秒

% 地面UE移动带来的多普勒
UE_max_speed_m_s = UE_max_speed*(1000/3600);% 单位m/s
Doppler_UE_max = (UE_max_speed_m_s*car_frequency)/c; %单位Hz

% 时延抽头的尺度 % 单位 秒
one_delay_tap = 1/(M*delta_f);
% 多普勒抽头的尺度 % 单位 Hz
one_doppler_tap = 1/(N*delta_t);

%%%%%%%%%%%%%%%%%%%量化为数字信道（NTN-TDL-D）%%%%%%%%%%%%%%%%
% 多径参数，计划遵循NTN—TDL—D模型，（6.9.2，P91）
% 该模型为存在LOS地卫星模型，共有三条径
% 第一条为LOS+瑞利（即莱斯分布）（LOS概率见6.6.1节P48）
% 第二为低时延 径
% 第三为高时延 径


% 俯仰角为 50度 的数据采集

delay_normalized = [0 0 0.5596 7.3340]; % 3Gpp NTN-TDL-D 归一化的时延参数
pdp =[-0.284,-11.991,-9.887,-16.771];   % 城市信道的阴影衰落 power delay profile 功率-时延谱，单位dB（38.811 NTN—TDL—D）
delay_RMS = 250*10^-9;                  % 均方根时延扩展参数，单位秒



% 延迟的相关计算
T_rms = delay_RMS/one_delay_tap; % 归一化均方根时延扩展参数，用于LMMSE信道估计，不参与函数的计算
delay_taps = (delay_sat+delay_normalized*delay_RMS)/one_delay_tap; % 考虑了卫星的大时延
delay_taps = round(delay_taps);     % assuming no fraction for the delay，时延不计算分数
taps = length(delay_taps);          % number of delay taps，确认有多少个抽头
% 功率的相关计算
pow_prof_taps = 10.^(pdp/10);    % 将dB转化为正常功率增益
pow_prof_taps = pow_prof_taps/sum(pow_prof_taps); % 归一化能量，即总和相加为1
% 根据每条信道的信道参数，得到每条径的的幅度，并考虑信道的类型，即：3Gpp NTN-TDL-D模型的瑞利信道或AWGN信道
 chan_amplitude_coef_taps = sqrt(pow_prof_taps).*(sqrt(1/2) * (randn(1,taps)+1i*randn(1,taps)));% channel coef. for each path （coefficient系数），其余径为瑞利信道
 chan_amplitude_coef_taps(1) = sqrt(pow_prof_taps(1)); % 首径为AWGN信道
% chan_amplitude_coef_taps = sqrt(pow_prof_taps); % 将所有瑞利信道变成AWGN信道
% 多普勒的相关计算
max_UE_Doppler_tap = Doppler_UE_max/one_doppler_tap;
sat_Doppler_tap = doppler_sat/one_doppler_tap;
doppler_taps = sat_Doppler_tap+(max_UE_Doppler_tap*cos(2*pi*rand(1,taps)));% 考虑星地间大多普勒
% doppler_taps = round(doppler_taps);     % 假设多普勒也不计算分数
% Doppler_taps = (max_UE_Doppler_tap*cos(2*pi*rand(1,taps)));% 假设星地间大多普勒已经补偿

if strcmp(channel_module,'NTN_TDL')
    doppler_taps = sat_Doppler_tap*ones(1,taps);% 此时地面UE不存在移动
elseif strcmp(channel_module,'NTN_TDL_Moving')
    % 无额外改动
elseif strcmp(channel_module,'Multipath_delay')
    doppler_taps = zeros(1,taps); % 无多普勒
elseif strcmp(channel_module,'Singlpath_doppler')
    delay_taps = delay_taps(1); % 单径
    taps = length(delay_taps); 
    chan_amplitude_coef_taps = 1;
    doppler_taps = sat_Doppler_tap;
elseif strcmp(channel_module,'Multipath_doppler')
    delay_taps = zeros(1,taps); % 无时延
else
    msg = '信道环境错误，本函数仅适用于：NTN_TDL, NTN_TDL_Moving, Multipath_delay, Singlpath_doppler, Multipath_doppler';
    error(msg)
end

% 废弃文本：
% 我们假设在50度角的时候，通过表6.9.2-4反推卫星发射器增益（P48 LOS功率占比0.726）（表6.9.2-4）
% dis_sat_UE_50= sqrt(R_earth^2*sin(5/18*pi)^2+sat_alt^2+2*sat_alt*R_earth)-R_earth*sin(5/18*pi);
% 
% %卫星发射端增益
% sat_Tx_Pow = 32.45+20*log10(car_frequency/10^9)+20*log10(dis_sat_UE_50)-10*log(0.726)+UE_Rx_Noise;
% % %卫星LOS径的增益（俯仰角相关）
% % sat_LOS = sat_Tx_Pow+10*log(table(table(:,1) == round(ele_ang),5))-PL_LOS;
% % %卫星NLOS径的增益（俯仰角相关）
% % sat_NLOS = sat_Tx_Pow+10*log(table(table(:,1) == round(ele_ang),5))-PL_NLOS;
% %多径的功率-时延谱重新改写，保证LOS的主径能量与3Gpp 38.811的6.6.1节相对应
% pow_prof_LOS_perc = table(table(:,1) == round(sat_ele_ang),5);
% 
% pow_prof = 10.^(pdp/10);
% pow_prof(1) = pow_prof_LOS_perc*(sum(pow_prof)-pow_prof(1))/(1-pow_prof_LOS_perc);
% pow_prof = pow_prof/sum(pow_prof);%normalization of power delay profile %%%归一化
% % pdp = 10*log10(pow_prof);
% % pdp = [SF_LOS,SF_NLOS]+pdp;
% % pow_prof = 10.^(pdp/10);
% % pow_prof = pow_prof/sum(pow_prof);%%%%加入阴影后的归一化（有问题）
% chan_coef = sqrt(pow_prof).*(sqrt(1/2) * (randn(1,taps)+1i*randn(1,taps)));%channel coef. for each path （coefficient系数）
% chan_coef(1) = pow_prof(1);
%%用户在接收端的功率（50俯仰角时增益为0）
% UE_Rx_Pow = sat_Tx_Pow-FSPL-UE_Rx_Noise+10*log(0.726);
% delay_taps(:,:)=0;
% Doppler_taps(:,:) =0;
% chan_coef = [1,0,0,0];

end