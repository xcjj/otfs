clc
fprintf('%s%s',interval_mode,'-OTFS(N,M,QAM size)=')
display([N,M,QAM_mode]);
% display(num2str(serv_t_percent(1)),'service time');
display(i_frame,'Number of frames');
display(SNR_dB,'SNR');
display(simulation_frame,'simulation frame');
if exist('info_tx_OFDM_DT','var')
display(result_avg_ber_OFDM(1,:),'Average coded BER OFDM');
display(result_avg_fer_OFDM(1,:),'Average coded FER OFDM');
end
% % display(avg_no_of_iterations_MRC_OFDM,'Average number of iterations for the MRC turbo decoder(OFDM)');
if exist('info_tx_OTFS_DT','var')
display(result_avg_ber_OTFS(1,:),'Average coded BER OTFS');
display(result_avg_fer_OTFS(1,:),'Average coded FER OTFS');
end
% display(avg_no_of_iterations_MRC,'Average number of iterations for the MRC turbo decoder(OTFS)');


