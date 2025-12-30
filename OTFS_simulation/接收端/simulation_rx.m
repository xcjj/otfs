%% 仿真末端（接收部分）

% 高频解调（无）
% 低频解调（无）
% 信号同步（小时延） 上文在大时延处 已经进行了相应的简化，此处也不用编写 


% 去掉保护间隔/循环前缀，同时串并转换，这里默认知道矩阵尺寸, 输出尺寸（M，N）
 
    
    if exist('info_tx_OTFS_DT','var')

        if strcmp(interval_mode,'None')
            info_rx_OTFS_S2P_DT  = reshape(info_channel_out_OTFS_DT,size(info_tx_OTFS_interval_DT));
            info_rx_OTFS_deinterval_DT = info_rx_OTFS_S2P_DT;
        elseif strcmp(interval_mode,'CP')
            info_rx_OTFS_S2P_DT  = reshape(info_channel_out_OTFS_DT,size(info_tx_OTFS_interval_DT));
            info_rx_OTFS_deinterval_DT = info_rx_OTFS_S2P_DT((length_CP+1):size(info_tx_OTFS_interval_DT,1) , :);
        elseif strcmp(interval_mode,'ZP') 
            info_rx_OTFS_S2P_DT  = reshape(info_channel_out_OTFS_DT,size(info_tx_OTFS_interval_DT));
            info_rx_OTFS_deinterval_DT = info_rx_OTFS_S2P_DT(1:(size(info_tx_OTFS_interval_DT,1)-length_ZP) , :);

        elseif strcmp(interval_mode,'RCP')
            info_rx_OTFS_S2P_DT  = reshape(info_channel_out_OTFS_DT(length_RCP+1:size(info_channel_out_OTFS_DT,1)),size(info_tx_OTFS_interval_DT));
            info_rx_OTFS_deinterval_DT = info_rx_OTFS_S2P_DT;
        elseif strcmp(interval_mode,'RZP')
            info_rx_OTFS_S2P_DT  = reshape(info_channel_out_OTFS_DT(1:(size(info_channel_out_OTFS_DT,1)-length_RZP)),size(info_tx_OTFS_interval_DT));
            info_rx_OTFS_deinterval_DT = info_rx_OTFS_S2P_DT;
        end
    end



     if exist('info_tx_OFDM_DT','var')

        if strcmp(interval_mode,'None')
            info_rx_OFDM_S2P_DT  = reshape(info_channel_out_OFDM_DT,size(info_tx_OFDM_interval_DT));
            info_rx_OFDM_deinterval_DT = info_rx_OFDM_S2P_DT;
        elseif strcmp(interval_mode,'CP')
            info_rx_OFDM_S2P_DT  = reshape(info_channel_out_OFDM_DT,size(info_tx_OFDM_interval_DT));
            info_rx_OFDM_deinterval_DT = info_rx_OFDM_S2P_DT((length_CP+1):(M+length_CP) , :);
        elseif strcmp(interval_mode,'ZP') 
            info_rx_OFDM_S2P_DT  = reshape(info_channel_out_OFDM_DT,size(info_tx_OFDM_interval_DT));
            info_rx_OFDM_deinterval_DT = info_rx_OFDM_S2P_DT(1:M , :);

        elseif strcmp(interval_mode,'RCP')
            info_rx_OFDM_S2P_DT  = reshape(info_channel_out_OFDM_DT(length_RCP+1:size(info_channel_out_OFDM_DT,1)),size(info_tx_OFDM_interval_DT));
            info_rx_OFDM_deinterval_DT = info_rx_OFDM_S2P_DT;
        elseif strcmp(interval_mode,'RZP')
            info_rx_OFDM_S2P_DT  = reshape(info_channel_out_OFDM_DT(1:M*N),size(info_tx_OFDM_interval_DT));
            info_rx_OFDM_deinterval_DT = info_rx_OFDM_S2P_DT;
        end
    end

 




% 信道估计

if exist('info_tx_OTFS_DT','var')
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
end


if exist('info_tx_OFDM_DT','var')
    if strcmp(channel_estimation_mode_OFDM,'None')
        info_rx_OFDM_estimation_DT = info_rx_OFDM_deinterval_DT; 
    else
        [info_rx_OFDM_domain,channel_H_est_OFDM_domain] = ...
            function_channel_estimation...
            (M_pilot_OFDM,N_pilot_OFDM,Fn_M_pilot,Fn_N_pilot,...
            info_rx_OFDM_deinterval_DT,pilot_symbols_OFDM,...
            QAM_mode,pilot_form_OFDM,pilot_domain_OFDM,channel_estimation_mode_OFDM,...
            pilot_Indx_OFDM,SNR(i_SNR),sigma_2(i_SNR), Channel_HDT,interval_mode,length_interval);

    end
end

% 信道均衡,输出尺寸（M_data，N_data）

if exist('info_tx_OTFS_DT','var')
    if strcmp(channel_equalization_mode_OTFS,'None')
        info_rx_OTFS_equalization_DT = info_rx_OTFS_estimation_DT; 
    elseif strcmp(channel_equalization_mode_OTFS,'ZF')||...
            strcmp(channel_equalization_mode_OTFS,'LS')||...
            strcmp(channel_equalization_mode_OTFS,'MMSE')

        info_rx_OTFS_equalization_DT = function_channel_equalization...
            (Fn_M_OTFS,Fn_N_OTFS,Fn_M_pilot,Fn_N_pilot,...
            channel_equalization_mode_OTFS,pilot_domain_OTFS,pilot_form_OTFS,...
            info_rx_OTFS_domain,channel_H_est_OTFS_domain,...
            Data_Indx_OTFS,sigma_2(i_SNR));

    end
end

if exist('info_tx_OFDM_DT','var')
    if strcmp(channel_equalization_mode_OFDM,'None')
        info_rx_OFDM_equalization_DT = info_rx_OFDM_estimation_DT; 
    elseif strcmp(channel_equalization_mode_OFDM,'ZF')||...
            strcmp(channel_equalization_mode_OFDM,'LS')||...
            strcmp(channel_equalization_mode_OFDM,'MMSE')

        info_rx_OFDM_equalization_DT = function_channel_equalization...
            (Fn_M_OFDM,Fn_N_OFDM,Fn_M_pilot,Fn_N_pilot,...
            channel_equalization_mode_OFDM,pilot_domain_OFDM,pilot_form_OFDM,...
            info_rx_OFDM_domain,channel_H_est_OFDM_domain,...
            Data_Indx_OFDM,sigma_2(i_SNR));
    end
end


% 数字解调（OTFS/OFDM）,输出尺寸（M，N_data）

    if exist('info_tx_OTFS_DT','var') 
        info_rx_deOTFS_DD = info_rx_OTFS_equalization_DT * Fn_N_OTFS; % 归一化的FFT矩阵是酉矩阵，共轭转秩等于矩阵的逆
    end
    if exist('info_tx_OFDM_DT','var')
        info_rx_deOFDM_FT = Fn_M_OFDM * info_rx_OFDM_equalization_DT; 
    end
% 并串转换    输出尺寸（1,M*N_data）
    if exist('info_tx_OTFS_DT','var') 
        info_rx_OTFS_P2S = reshape(info_rx_deOTFS_DD,1,numel(info_rx_deOTFS_DD));
    end
    if exist('info_tx_OFDM_DT','var')
        info_rx_OFDM_P2S = reshape(info_rx_deOFDM_FT,1,numel(info_rx_deOFDM_FT)); 
    end


% 解QAM，输出尺寸（QAM*M*N_data,1）
if strcmp(channel_coding_mode,'None')
    if exist('info_tx_OTFS_DT','var') 
        info_rx_bits_OTFS_deQAM = qamdemod(info_rx_OTFS_P2S,QAM_mode,'gray','OutputType','bit');
        info_rx_bits_OTFS_deQAM = reshape(info_rx_bits_OTFS_deQAM,numel(info_rx_bits_OTFS_deQAM),1);
    end
    if exist('info_tx_OFDM_DT','var')
        info_rx_bits_OFDM_deQAM = qamdemod(info_rx_OFDM_P2S,QAM_mode,'gray','OutputType','bit');
        info_rx_bits_OFDM_deQAM = reshape(info_rx_bits_OFDM_deQAM,numel(info_rx_bits_OFDM_deQAM),1);
    end
else
    if exist('info_tx_OTFS_DT','var') 
        info_rx_bits_OTFS_deQAM = qamdemod(info_rx_OTFS_P2S,QAM_mode,'gray','OutputType','llr');
        info_rx_bits_OTFS_deQAM = reshape(info_rx_bits_OTFS_deQAM,numel(info_rx_bits_OTFS_deQAM),1);
    end
    if exist('info_tx_OFDM_DT','var')
        info_rx_bits_OFDM_deQAM = qamdemod(info_rx_OFDM_P2S,QAM_mode,'gray','OutputType','llr');
        info_rx_bits_OFDM_deQAM = reshape(info_rx_bits_OFDM_deQAM,numel(info_rx_bits_OFDM_deQAM),1);
    end
end

% 解交织，输出尺寸（QAM*M*N_data,1）
    if strcmp(interleaving_mode, 'Interleaving')
        if exist('info_tx_OTFS_DT','var') 
            info_rx_bits_OTFS_deinterleaved = info_rx_bits_OTFS_deQAM(RND_Intrlv_rev); 
        end
        if exist('info_tx_OFDM_DT','var')
            info_rx_bits_OFDM_deinterleaved = info_rx_bits_OFDM_deQAM(RND_Intrlv_rev); 
        end
    else
        if exist('info_tx_OTFS_DT','var') 
            info_rx_bits_OTFS_deinterleaved = info_rx_bits_OTFS_deQAM; 
        end
        if exist('info_tx_OFDM_DT','var')
            info_rx_bits_OFDM_deinterleaved = info_rx_bits_OFDM_deQAM; 
        end
    end

% 解码
    if exist('info_tx_OTFS_DT','var') 
        if strcmp(channel_coding_mode,'None')
            info_rx_bits_OTFS_decoded = info_rx_bits_OTFS_deinterleaved;    % 尺寸为MxNxQAM
        elseif strcmp(channel_coding_mode,'LDPC')
            info_rx_bits_OTFS_decoded = zeros(info_bits_length,1);% 初始化数据帧，没有码元的部分被置零，尺寸为MxNxQAM
            for i_LDPC_block=1:1:LDPC_trans_blocks
                t1 = (i_LDPC_block-1)*LDPC_info_length + 1;     % 确定编码后码块在数据帧中的位置
                t2 = (i_LDPC_block-1)*LDPC_codeword_length + 1; % 确定数据码元在数据帧中的位置
                info_rx_bits_OTFS_decoded(t1:t1+LDPC_info_length-1,1) = step(hDec, info_rx_bits_OTFS_deinterleaved(t2:t2+LDPC_codeword_length-1,1)); %解码
            end
        elseif strcmp(channel_coding_mode,'Polar')
            info_rx_bits_OTFS_decoded = zeros(info_bits_length,1);% 初始化数据帧，没有码元的部分被置零，尺寸为MxNxQAM
            for i_Polar_block=1:1:Polar_trans_blocks
                t1 = (i_Polar_block-1)*Polar_info_length + 1;     % 确定编码后码块在数据帧中的位置
                t2 = (i_Polar_block-1)*Polar_codeword_length + 1; % 确定数据码元在数据帧中的位置
                info_rx_bits_OTFS_decoded(t1:t1+Polar_info_length-1,1) = nrPolarDecode(info_rx_bits_OTFS_deinterleaved(t2:t2+Polar_codeword_length-1,1),Polar_info_length,Polar_codeword_length,Polar_CRC_decode_length);%解码
            end
        end
    end
    if exist('info_tx_OFDM_DT','var') 
        if strcmp(channel_coding_mode,'None')
            info_rx_bits_OFDM_decoded = info_rx_bits_OFDM_deinterleaved;    % 尺寸为MxNxQAM
        elseif strcmp(channel_coding_mode,'LDPC')
            info_rx_bits_OFDM_decoded = zeros(info_bits_length,1);% 初始化数据帧，没有码元的部分被置零，尺寸为MxNxQAM
            for i_LDPC_block=1:1:LDPC_trans_blocks
                t1 = (i_LDPC_block-1)*LDPC_info_length + 1;     % 确定编码后码块在数据帧中的位置
                t2 = (i_LDPC_block-1)*LDPC_codeword_length + 1; % 确定数据码元在数据帧中的位置
                info_rx_bits_OFDM_decoded(t1:t1+LDPC_info_length-1,1) = step(hDec, info_rx_bits_OFDM_deinterleaved(t2:t2+LDPC_codeword_length-1,1)); %解码
            end
        elseif strcmp(channel_coding_mode,'Polar')
             info_rx_bits_OFDMs_decoded = zeros(info_bits_length,1);% 初始化数据帧，没有码元的部分被置零，尺寸为MxNxQAM
            for i_Polar_block=1:1:Polar_trans_blocks
                t1 = (i_Polar_block-1)*Polar_info_length + 1;     % 确定编码后码块在数据帧中的位置
                t2 = (i_Polar_block-1)*Polar_codeword_length + 1; % 确定数据码元在数据帧中的位置
                info_rx_bits_OFDM_decoded(t1:t1+Polar_info_length-1,1) = nrPolarDecode(info_rx_bits_OFDM_deinterleaved(t2:t2+Polar_codeword_length-1,1),Polar_info_length,Polar_codeword_length,Polar_CRC_decode_length);%解码
            end
        end
    end
