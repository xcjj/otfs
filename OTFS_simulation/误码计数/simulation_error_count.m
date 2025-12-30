
% 误码计数
    if exist('info_tx_OFDM_DT','var') 
        %% errors count%%%%%   (OFDM)  
        errors_OFDM = sum(xor(info_rx_bits_OFDM_decoded,info_tx_bits));
        if(errors_OFDM>0)
            result_err_ber_OFDM(1,i_SNR) = result_err_ber_OFDM(1,i_SNR) + errors_OFDM;           % Bit error rate
            result_err_fer_OFDM(1,i_SNR) = result_err_fer_OFDM(1,i_SNR) + 1;                    % Frame error rate: frame error happends if even one bit in a frame is wronf
        end
        result_avg_ber_OFDM(1,i_SNR) = result_err_ber_OFDM(1,i_SNR).'/length(info_tx_bits)/i_frame;
        result_avg_fer_OFDM(1,i_SNR) = result_err_fer_OFDM(1,i_SNR).'/i_frame;

    end
        %% errors count%%%%%   (OTFS)  
    if exist('info_tx_OTFS_DT','var') 
        errors_OTFS = sum(xor(info_rx_bits_OTFS_decoded,info_tx_bits));
        if(errors_OTFS>0)
            result_err_ber_OTFS(1,i_SNR) = result_err_ber_OTFS(1,i_SNR) + errors_OTFS;           % Bit error rate
            result_err_fer_OTFS(1,i_SNR) = result_err_fer_OTFS(1,i_SNR) + 1;                    % Frame error rate: frame error happends if even one bit in a frame is wronf
        end
        result_avg_ber_OTFS(1,i_SNR) = result_err_ber_OTFS(1,i_SNR).'/length(info_tx_bits)/i_frame;
        result_avg_fer_OTFS(1,i_SNR) = result_err_fer_OTFS(1,i_SNR).'/i_frame;

    end

    simulation_frame(1,i_SNR) = i_frame;
