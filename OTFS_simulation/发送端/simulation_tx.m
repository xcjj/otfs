%% 仿真前端（发送端部分）

% 信号生成，列向量，长度不固定
    info_tx_bits = randi([0,1],info_bits_length,1); % 数据长度取决于所用编码，以及导频
    % info_tx_bits = zeros(info_bits_length,1); % 数据长度取决于所用编码，以及导频
% 编码，输出尺寸（M*N_data*QAM，1）
    if strcmp(channel_coding_mode,'None')
        info_tx_bits_coded = info_tx_bits;    % 尺寸为M_data*N_data*QAM
    elseif strcmp(channel_coding_mode,'LDPC')
        info_tx_bits_coded = zeros(bits_data_perframe,1);% 初始化数据帧，没有码元的部分被置零，尺寸为M_data*N_data*QAM
        for i_LDPC_block=1:1:LDPC_trans_blocks
            t1 = (i_LDPC_block-1)*LDPC_info_length + 1;     % 确定编码后码块在数据帧中的位置
            t2 = (i_LDPC_block-1)*LDPC_codeword_length + 1; % 确定数据码元在数据帧中的位置
            info_tx_bits_coded(t2:t2+(LDPC_codeword_length)-1,1) = step(hEnc, info_tx_bits(t1:t1+LDPC_info_length-1));% 将LDPC码块嵌入数据帧中
        end
    elseif strcmp(channel_coding_mode,'Polar')
        info_tx_bits_coded = zeros(bits_data_perframe,1);% 初始化数据帧，没有码元的部分被置零，尺寸为M_data*N_data*QAM
        for i_Polar_block=1:1:Polar_trans_blocks
            t1 = (i_Polar_block-1)*Polar_info_length + 1;     % 确定编码后码块在数据帧中的位置
            t2 = (i_Polar_block-1)*Polar_codeword_length + 1; % 确定数据码元在数据帧中的位置
            info_tx_bits_coded(t2:t2+(Polar_codeword_length)-1,1) = nrPolarEncode(info_tx_bits(t1:t1+Polar_info_length-1),Polar_codeword_length);% 将Polar码块嵌入数据帧中
        end
    end

% 交织(随机交织)，输出尺寸（M_data*N_data*QAM，1）
    if strcmp(interleaving_mode, 'Interleaving')
        state=randi(10000);
        RND_Intrlv     = randintrlv  (1:bits_data_perframe,state);   % 交织模块
        RND_Intrlv_rev = randdeintrlv(1:bits_data_perframe,state);   % 解交织模块
        info_tx_bits_interleaved = info_tx_bits_coded(RND_Intrlv); 
    elseif strcmp(interleaving_mode, 'None')
        info_tx_bits_interleaved = info_tx_bits_coded;
    else
        msg = '选择交织模式: (None/Interleaving)';
        error(msg)
    end

% 调幅调相（QAM）输出尺寸（1, M_data*N_data）
        info_tx_symbols_QAM=qammod(reshape(info_tx_bits_interleaved,QAM_bits,symbols_data_perframe), QAM_mode,'gray','InputType','bit');  %格雷码  
% 串并转换 sreializer to parallel ,输出尺寸 (M_data,N_data)
        info_tx_symbols_S2P_OTFS = reshape(info_tx_symbols_QAM,M_data_OTFS,N_data_OTFS); % 按列读取，按列摆放
        info_tx_symbols_S2P_OFDM = reshape(info_tx_symbols_QAM,M_data_OFDM,N_data_OFDM); % 按列读取，按列摆放
% 数字调制（OTFS/OFDM）,输出尺寸 (M_data,N_data)
    if strcmp(Modem_mode_1,'OTFS')||strcmp(Modem_mode_2,'OTFS')
        info_tx_OTFS_DT = info_tx_symbols_S2P_OTFS*Fn_N_OTFS';% DT代表delay_time域 
    end
    if strcmp(Modem_mode_1,'OFDM')||strcmp(Modem_mode_2,'OFDM')
        info_tx_OFDM_DT = Fn_M_OFDM'*info_tx_symbols_S2P_OFDM; 
    end


% 插入导频（频域导频），输出尺寸 (M,N)

    if exist('info_tx_OTFS_DT','var')
        pilot_symbols_OTFS = function_pilot_generation(QAM_mode,M_data_OTFS,N_data_OTFS,pilot_form_OTFS,embeded_pilot_amp);
        if strcmp(channel_estimation_mode_OTFS,'None')
            info_tx_OTFS_pilot_DT = info_tx_OTFS_DT;
        else
            info_tx_OTFS_pilot_DT = function_pilot_insert(M,N,Fn_M_OTFS,Fn_N_OTFS,Fn_M_pilot,Fn_N_pilot,info_tx_OTFS_DT,pilot_symbols_OTFS,pilot_domain_OTFS,pilot_form_OTFS,Data_Indx_OTFS,pilot_Indx_OTFS);
        end
    end
    if exist('info_tx_OFDM_DT','var')
        pilot_symbols_OFDM = function_pilot_generation(QAM_mode,M_data_OFDM,N_data_OFDM,pilot_form_OFDM,embeded_pilot_amp);
        if strcmp(channel_estimation_mode_OFDM,'None')
            info_tx_OFDM_pilot_DT = info_tx_OFDM_DT;
        else
            info_tx_OFDM_pilot_DT = function_pilot_insert(M,N,Fn_M_OFDM,Fn_N_OFDM,Fn_M_pilot,Fn_N_pilot,info_tx_OFDM_DT,pilot_symbols_OFDM,pilot_domain_OFDM,pilot_form_OFDM,Data_Indx_OFDM,pilot_Indx_OFDM);
        end

    end




% 插入保护间隔（循环前缀），CP加在前面，ZP加在后面

    if exist('info_tx_OTFS_DT','var')
        if strcmp(interval_mode,'None') % 输出尺寸（M，N）
            info_tx_OTFS_interval_DT = info_tx_OTFS_pilot_DT;
            info_tx_OTFS_RCP_DT  = [];
            info_tx_OTFS_RZP_DT  = [];
        elseif strcmp(interval_mode,'CP') % 输出尺寸（Length_CP+M，N）
            info_tx_OTFS_interval_DT = [info_tx_OTFS_pilot_DT((M-length_CP+1):M,1:N);info_tx_OTFS_pilot_DT];
            info_tx_OTFS_RCP_DT  = [];
            info_tx_OTFS_RZP_DT  = [];
        elseif strcmp(interval_mode,'ZP') % 输出尺寸 (M+Length_ZP，N)
            info_tx_OTFS_interval_DT = [info_tx_OTFS_pilot_DT;zeros(length_ZP,N)];
            info_tx_OTFS_RCP_DT  = [];
            info_tx_OTFS_RZP_DT  = [];
        elseif strcmp(interval_mode,'RCP') % 输出尺寸 (length_CP,1)+（M，N）
            info_tx_OTFS_interval_DT = info_tx_OTFS_pilot_DT;
            info_tx_OTFS_RCP_DT  = info_tx_OTFS_pilot_DT((M-length_RCP+1):M,N);
            info_tx_OTFS_RZP_DT  = [];
        elseif strcmp(interval_mode,'RZP') % 输出尺寸 (M，N)+(length_ZP,1)
            info_tx_OTFS_interval_DT = info_tx_OTFS_pilot_DT;
            info_tx_OTFS_RCP_DT  = [];
            info_tx_OTFS_RZP_DT  = info_tx_OTFS_pilot_DT((M-length_RZP+1):M,N);
        end

    end

    if exist('info_tx_OFDM_DT','var')
        if strcmp(interval_mode,'None') % 输出尺寸（M，N）
            info_tx_OFDM_interval_DT = info_tx_OFDM_pilot_DT;
            info_tx_OFDM_RCP_DT  = [];
            info_tx_OFDM_RZP_DT  = [];
        elseif strcmp(interval_mode,'CP') % 输出尺寸（Length_CP+M，N）
            info_tx_OFDM_interval_DT = [info_tx_OFDM_pilot_DT((M-length_CP+1):M,1:N);info_tx_OFDM_pilot_DT];
            info_tx_OFDM_RCP_DT  = [];
            info_tx_OFDM_RZP_DT  = [];
        elseif strcmp(interval_mode,'ZP') % 输出尺寸 (M+Length_ZP，N)
            info_tx_OFDM_interval_DT = [info_tx_OFDM_pilot_DT;zeros(length_ZP,N)];
            info_tx_OFDM_RCP_DT  = [];
            info_tx_OFDM_RZP_DT  = [];
        elseif strcmp(interval_mode,'RCP') % 输出尺寸 (length_CP,1)+（M，N）
            info_tx_OFDM_interval_DT = info_tx_OFDM_pilot_DT;
            info_tx_OFDM_RCP_DT  = info_tx_OFDM_pilot_DT((M-length_RCP+1):M,N);
            info_tx_OFDM_RZP_DT  = [];
        elseif strcmp(interval_mode,'RZP') % 输出尺寸 (M，N)+(length_ZP,1)
            info_tx_OFDM_interval_DT = info_tx_OFDM_pilot_DT;
            info_tx_OFDM_RCP_DT  = [];
            info_tx_OFDM_RZP_DT  = info_tx_OFDM_pilot_DT((M-length_RZP+1):M,N);
        end

    end


% 并串转换，parallel to sreializer ,输出尺寸 (MxN,1),或者((M+CP/ZP)xN,1),或者(MxN+RCP/RZP,1)
    
    if exist('info_tx_OTFS_DT','var')    
        info_tx_OTFS_P2S_DT = [info_tx_OTFS_RCP_DT; reshape(info_tx_OTFS_interval_DT,numel(info_tx_OTFS_interval_DT),1); info_tx_OTFS_RZP_DT]; 
    end
    if exist('info_tx_OFDM_DT','var')    
        info_tx_OFDM_P2S_DT = [info_tx_OFDM_RCP_DT; reshape(info_tx_OFDM_interval_DT,numel(info_tx_OFDM_interval_DT),1); info_tx_OFDM_RZP_DT]; 
    end

% 低频调制（无）
% 高频调制（无）
% MIMO（无）
