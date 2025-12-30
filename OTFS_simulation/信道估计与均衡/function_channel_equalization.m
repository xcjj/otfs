function info_rx_data_equalization_DT = function_channel_equalization...
(Fn_M,Fn_N,Fn_M_pilot,Fn_N_pilot,channel_equalization_mode,pilot_domain,pilot_form,info_rx_domain,channel_H_estimation,data_Indx,sigma_2)
%FUNCTION_CHANNEL_EQUALIZATION 此处显示有关此函数的摘要
%   此处显示详细说明
%   该均衡代码传入的信号为带导频的信号，信号进行均衡后再提取出数据信号

% 域的转换，取决于所选调制域，统一转换为FT域进行均衡
    if strcmp(pilot_domain, 'DT')
        channel_H_estimation_FT = Fn_M_pilot*channel_H_estimation; 
        info_rx_FT= Fn_M_pilot*info_rx_domain;
    elseif strcmp(pilot_domain, 'FT')
        channel_H_estimation_FT = channel_H_estimation; 
        info_rx_FT= info_rx_domain;
    elseif strcmp(pilot_domain, 'DD')
        if contains(channel_equalization_mode,'OTFS_Optimize')

        else
            channel_H_estimation_FT = Fn_M_pilot*channel_H_estimation*Fn_N_pilot'; 
            info_rx_FT= Fn_M_pilot*info_rx_domain*Fn_N_pilot';
        end
    elseif strcmp(pilot_domain, 'FD')
        channel_H_estimation_FT = channel_H_estimation*Fn_N_pilot'; 
        info_rx_FT= info_rx_domain*Fn_N_pilot';
    else 
        msg = '请选择正确的导频插入域 （DT，FT, DD, FD）,或检测变量“channel_equalization_mode”包含的“OTFS_Optimize”拼写是否错误';
        error(msg)
    end

% 带导频的均衡
    if strcmp(channel_equalization_mode,'LS')
        info_rx_equalization_FT=info_rx_FT.*conj(channel_H_estimation_FT)./(abs(channel_H_estimation_FT).^2);
    elseif strcmp(channel_equalization_mode,'ZF')
        info_rx_equalization_FT=info_rx_FT./conj(channel_H_estimation_FT);       
    elseif strcmp(channel_equalization_mode,'MMSE')
        % info_rx_equalization_domain=info_rx_data_domain.*conj(channel_H_data_estimation)./(abs(channel_H_data_estimation).^2);
        msg = 'MMSE均衡没做好，选择LS均衡即可';
        error(msg)
    else
        msg = '暂不支持其它形式，请选择合适的信道均衡方式（LS，ZF,MMSE）';
        error(msg)
    end

% 转换回对应的域，并提取出数据

    if strcmp(pilot_domain, 'DT')
        info_rx_equalization_domain = Fn_M_pilot'*info_rx_equalization_FT;
    elseif strcmp(pilot_domain, 'FT')
        info_rx_equalization_domain = info_rx_equalization_FT;
    elseif strcmp(pilot_domain, 'DD')
        info_rx_equalization_domain = Fn_M_pilot'*info_rx_equalization_FT*Fn_N_pilot;
    elseif strcmp(pilot_domain, 'FD')
        info_rx_equalization_domain = info_rx_equalization_FT*Fn_N_pilot';
    end


    if strcmp(pilot_form, 'block')
        info_rx_data_equalization_domain = info_rx_equalization_domain(:,data_Indx);
    elseif strcmp(pilot_form, 'comb_like')
        info_rx_data_equalization_domain = info_rx_equalization_domain(data_Indx,:);
    elseif strcmp(pilot_form, 'lattice')
        msg = '暂不支持格状导频';
        error(msg)
    elseif strcmp(pilot_form, 'OTFS_delta')
        info_rx_data_equalization_domain = info_rx_equalization_domain(data_Indx,:);
    end


% 最后将信号转换成DT信号
    if strcmp(pilot_domain, 'DT')
        info_rx_data_equalization_DT = info_rx_data_equalization_domain;
    elseif strcmp(pilot_domain, 'FT')
        info_rx_data_equalization_DT = Fn_M'*info_rx_data_equalization_domain; 
    elseif strcmp(pilot_domain, 'DD')
        info_rx_data_equalization_DT = info_rx_data_equalization_domain*Fn_N'; 
    elseif strcmp(pilot_domain, 'FD')
        info_rx_data_equalization_DT = Fn_M'*info_rx_data_equalization_domain*Fn_N';    
    
    end



end

