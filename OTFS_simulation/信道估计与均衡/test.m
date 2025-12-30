



% 部分参数初始化
[M,N] = size(info_rx_OTFS_deinterval_DT);

if strcmp(channel_estimation_mode,'LMMSE')
    % LMMSE估计，此估计无法用于估计单径信道，
    % trms为多经信道的平均延时，t_max为最大延时,此处所有的时间都是已经对采样间隔做了归一化后的结果
    % beta = E(x^2)*E(1/x^2); x代表QAM信号的星座值 https://blog.csdn.net/shenyuhou/article/details/121871048
    if QAM_mode == 64
        beta=2.6857;
    elseif  QAM_mode == 16
        beta=17/9;
    elseif  QAM_mode == 4
        beta=1;
    else
        msg = '本代码的LMMSE需要QAM值在4，16，64中选择，其余未计算，如有需要，找到对应部分更改相关代码';
        error(msg)
    end

%     if sum(strcmp(channel_module,{'NTN_TDL','NTN_TDL_Moving','Multipath_delay'}))
%     trms=T_rms;
%     t_max=max(channel_delay_taps);
%     else
%     trms=1e-6;
%     t_max=1e-6; 
%     end
%     Rhh=zeros(M,M); %协方差矩阵,利用统计信息生成信道的自相关矩阵
%     for k=1:M
%         for l=1:M
%             Rhh(k,l)=(1-exp((-1)*t_max*((1/trms)+1i*2*pi*(k-l)/M)))./(trms*(1-exp((-1)*t_max/trms))*((1/trms)+1i*2*pi*(k-l)/M));
%         end
%     end
end


% 域的转换，取决于所选调制域
if strcmp(pilot_domain_OTFS, 'DT')
    % 不需要更改
    info_rx_domain = info_rx_OTFS_deinterval_DT;
elseif strcmp(pilot_domain_OTFS, 'FT')
    info_rx_domain = Fn_M_pilot*info_rx_OTFS_deinterval_DT; % 信号转换为FT域
elseif strcmp(pilot_domain_OTFS, 'DD')
    info_rx_domain = info_rx_OTFS_deinterval_DT*Fn_N_pilot; % 信号转换为DD域
else 
    msg = '请选择正确的导频插入域 （DT，FT, DD）';
    error(msg)
end


if strcmp(pilot_form, 'block')
    % 接收信号的导频数据与码元提取
    info_rx_pilot_domain = info_rx_domain(:,pilot_Indx);
    info_rx_data_domain  = info_rx_domain(:,Data_Indx);
    % 导频图案
    pilot_pattren = repmat(pilot_symbols,1,N_pilot); 
    if (strcmp(channel_estimation_mode,'LS')||strcmp(channel_estimation_mode,'LMMSE'))
        % 导频位置的信道响应LS估计,LS,LMMSE都要用
            pilot_estimation_domain=info_rx_pilot_domain./pilot_pattren; % LS 导频估计
        if strcmp(channel_estimation_mode,'LMMSE')
            Rhh = pilot_estimation_domain*pilot_estimation_domain';% 协方差矩阵
            Q=Rhh*(eye(M)/(Rhh+(beta/SNR(i_SNR))*eye(M))); %eye()斜对角为1
            pilot_estimation_domain=Q*pilot_estimation_domain; % LMMSE的估计矩阵
        end
        % 线性插值
        channel_H_domain = zeros(M,N);
        for i=1:M
            channel_H_domain(i,:)=interp1(pilot_Indx,pilot_estimation_domain(i,1:(N_pilot)),1:N,'linear');
        end
        channel_H_data_estimation=channel_H_domain(:,Data_Indx);
    end
    

elseif strcmp(pilot_form, 'lattice')
    info_rx_pilot_domain = info_rx_domain(pilot_Indx,:);
    info_rx_data_domain  = info_rx_domain(Data_Indx,:);
    pilot_pattren = repmat(pilot_symbols,M_pilot,1); % 导频图案

    if (strcmp(channel_estimation_mode,'LS')||strcmp(channel_estimation_mode,'LMMSE'))
        % 导频位置的信道响应LS估计,LS,LMMSE都要用
            pilot_estimation_domain=info_rx_pilot_domain./pilot_pattren; % LS 导频估计
        if strcmp(channel_estimation_mode,'LMMSE')
            Rhh = pilot_estimation_domain'*pilot_estimation_domain;% 协方差矩阵
            Q=Rhh*(eye(N)/(Rhh+(beta/SNR(i_SNR))*eye(N))); %eye()斜对角为1
            pilot_estimation_domain=pilot_estimation_domain*Q'; % LMMSE的估计矩阵
        end
        % 线性插值
        channel_H_domain = zeros(M,N);
        for j=1:N
            channel_H_domain(:,j)=interp1(pilot_Indx,pilot_estimation_domain(1:(M_pilot)),1:M,'linear');
        end
        channel_H_data_estimation=channel_H_domain(Data_Indx,:);
    end



elseif strcmp(pilot_form, 'comb_like')
    msg = '暂不支持梳状导频';
    error(msg)
end

