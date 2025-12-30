function [info_rx_domain,channel_H_domain] = function_channel_estimation ...
(M_pilot,N_pilot,Fn_M_pilot,Fn_N_pilot,info_rx_DT,pilot_symbols,QAM_mode,pilot_form,pilot_domain,...
channel_estimation_mode,pilot_Indx,SNR,sigma_2,Channel_HDT,interval_mode,length_interval)
%FUNCTION_CHANNEL_ESTIMATION 此处显示有关此函数的摘要
%   此处显示详细说明



% 部分参数初始化
[M,N] = size(info_rx_DT);
M_data = M-M_pilot;
N_data = N-N_pilot;

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
if strcmp(pilot_domain, 'DT')
    % 不需要更改
    info_rx_domain = info_rx_DT;
elseif strcmp(pilot_domain, 'FT')
    info_rx_domain = Fn_M_pilot*info_rx_DT; % 信号转换为FT域
elseif strcmp(pilot_domain, 'DD')
    info_rx_domain = info_rx_DT*Fn_N_pilot; % 信号转换为DD域
    
elseif strcmp(pilot_domain, 'FD')
    info_rx_domain = Fn_M_pilot*info_rx_DT*Fn_N_pilot; % 信号转换为FD域    
else 
    msg = '请选择正确的导频插入域 （DT，FT, DD，FD）';
    error(msg)
end


if strcmp(channel_estimation_mode,'ideal')


%     Channel_HDT = reshape(Channel_HDT,size(Channel_DT,2)/N,N);
     Channel_HDT = Channel_HDT*sqrt(M);
    if strcmp(interval_mode,'None')
        channel_H_dt = Channel_HDT;
    elseif strcmp(interval_mode,'CP')
%         channel_H_dt = Channel_HDT((length_interval+1):size(Channel_HDT,1) , :);
        channel_H_dt = Channel_HDT(1:(size(Channel_HDT,1)-length_interval) , :);
    elseif strcmp(interval_mode,'ZP') 
        channel_H_dt = Channel_HDT(1:(size(Channel_HDT,1)-length_interval) , :);
    elseif strcmp(interval_mode,'RCP')
        channel_H_dt = Channel_HDT;
    elseif strcmp(interval_mode,'RZP')
        channel_H_dt = Channel_HDT;
    end
        
    if strcmp(pilot_domain, 'DT')
        % 不需要更改
        channel_H_domain = channel_H_dt;
    elseif strcmp(pilot_domain, 'FT')
        channel_H_domain = Fn_M_pilot*channel_H_dt; % 信号转换为FT域
    elseif strcmp(pilot_domain, 'DD')
        channel_H_domain = channel_H_dt*Fn_N_pilot; % 信号转换为DD域
    elseif strcmp(pilot_domain, 'FD')
        channel_H_domain = Fn_M_pilot*channel_H_dt*Fn_N_pilot; % 信号转换为FD域    
    else 
        msg = '请选择正确的导频插入域 （DT，FT, DD，FD）';
        error(msg)
    end


elseif strcmp(pilot_form, 'block')
    % 接收信号的导频数据与码元提取
    info_rx_pilot_domain = info_rx_domain(:,pilot_Indx);
    % info_rx_data_domain  = info_rx_domain(:,data_Indx);
    % 导频图案
    pilot_pattren = repmat(pilot_symbols,1,N_pilot); 
    if (strcmp(channel_estimation_mode,'LS')||strcmp(channel_estimation_mode,'LMMSE'))
        % 导频位置的信道响应LS估计,LS,LMMSE都要用
            pilot_estimation_domain=info_rx_pilot_domain./pilot_pattren; % LS 导频估计
        if strcmp(channel_estimation_mode,'LMMSE')
            Rhh = pilot_estimation_domain*pilot_estimation_domain';% 协方差矩阵
            Q=Rhh*(eye(M)/(Rhh+(beta/SNR)*eye(M))); %eye()斜对角为1
            pilot_estimation_domain=Q*pilot_estimation_domain; % LMMSE的估计矩阵
        end
        % 线性插值
        channel_H_domain = zeros(M,N);
        for i=1:M
            channel_H_domain(i,:)=interp1(pilot_Indx,pilot_estimation_domain(i,1:(N_pilot)),1:N,'linear');
        end
        % channel_H_data_estimation_domain=channel_H_domain(:,data_Indx);
    end
    

elseif strcmp(pilot_form, 'comb_like')
    info_rx_pilot_domain = info_rx_domain(pilot_Indx,:);
    % info_rx_data_domain  = info_rx_domain(data_Indx,:);
    pilot_pattren = repmat(pilot_symbols,M_pilot,1); % 导频图案

    if (strcmp(channel_estimation_mode,'LS')||strcmp(channel_estimation_mode,'LMMSE'))
        % 导频位置的信道响应LS估计,LS,LMMSE都要用
            pilot_estimation_domain=info_rx_pilot_domain./pilot_pattren; % LS 导频估计
        if strcmp(channel_estimation_mode,'LMMSE')
            Rhh = pilot_estimation_domain'*pilot_estimation_domain;% 协方差矩阵
            Q=Rhh*(eye(N)/(Rhh+(beta/SNR)*eye(N))); %eye()斜对角为1
            pilot_estimation_domain=pilot_estimation_domain*Q'; % LMMSE的估计矩阵
        end
        % 线性插值
        channel_H_domain = zeros(M,N);
        for j=1:N
            channel_H_domain(:,j)=interp1(pilot_Indx,pilot_estimation_domain(1:(M_pilot),j),1:M,'linear');
        end
        % channel_H_data_estimation_domain=channel_H_domain(data_Indx,:);
    end



elseif strcmp(pilot_form, 'lattice')
    msg = '暂不支持格状导频';
    error(msg)


    
elseif strcmp(pilot_form, 'OTFS_delta')
    if strcmp(channel_estimation_mode,'embedded_pilot')&&strcmp(pilot_domain, 'DD')
        % 接收信号的导频数据与码元提取
        info_rx_pilot_domain = info_rx_domain(pilot_Indx,:);
        % info_rx_data_domain  = info_rx_domain(data_Indx,:);
        info_rx_pilot_domain(abs(info_rx_pilot_domain)<3*sqrt(sigma_2)) = 0; % 阈值判决

        pilot_estimation_domain = info_rx_pilot_domain(ceil(M_pilot/2):M_pilot,:)/pilot_symbols(ceil(N_data/2));% 接收端导频除以发射端导频能量，能量归一化
        pilot_estimation_domain = [pilot_estimation_domain(:,ceil(N_data/2):N_data),pilot_estimation_domain(:,1:(ceil(N_data/2)-1))];% 循环位移

        channel_H_domain = zeros(M,N);
        channel_H_domain(1:size(pilot_estimation_domain,1),1:size(pilot_estimation_domain,2)) = pilot_estimation_domain*sqrt(M*N);
        % channel_H_domain(abs(channel_H_domain*pilot_symbols(ceil(N_data/2))/N)<sqrt(sigma_2)) = 0;
        % channel_H_domain(abs(channel_H_domain)<sqrt(sigma_2)) = 0;
        % channel_H_domain(abs(channel_H_domain)<(max(max(abs(channel_H_domain)))/10+sqrt(sigma_2)/2)) = 0;


        % channel_H_data_estimation_domain=channel_H_domain(data_Indx,:);


        if(channel_H_domain == zeros(M,N))
            % 低信噪比情况下，有效符号可能会低于阈值，进而导致信道估计报错，为了防止报错，当出现这种情况时，将估计到的信道变成恒定矩阵
            channel_H_domain(1,1) = N;
        end
    else
        msg = '信道估计处使用embedded_pilot，但导频不是OTFS_delata格式，或不在DD域';
        error(msg)
    end

    
end



end

