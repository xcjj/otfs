% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 图表绘制 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 仿真完成，误码性能展示



dot =['x';'o';'s';'d';'h';'.'];
% color = ['r';'k';'b';'c'];

figure
    for i= 1:6
    eval(['semilogy(SNR_dB,ber_count_OTFS(i,:),''-',(dot(mod(i,6)+1,:)),''',''LineWidth'',2,''MarkerSize'',8)' ])
%      legend_str2{i}=['StudyCase ',Study_case,' with ',doppler_mode,' doppler in ',num2str(car_frequency/10^9),'GHz'];
    hold on
    end
%     for i= 1:3
%     eval(['semilogy(SNR_dB,ber_count_OFDM(i,:),''--',color(mod(i,4)+1,:),(dot(mod(i,4)+1,:)),''',''LineWidth'',2,''MarkerSize'',8)' ])
% %      legend_str2{i}=['StudyCase ',Study_case,' with ',doppler_mode,' doppler in ',num2str(car_frequency/10^9),'GHz'];
%     hold on
%     end
%     legend(legend_str2(9:16))
    axis([SNR_dB(1) SNR_dB(i_SNR) 1e-4 1e-0])
    grid on
    xlabel('SNR(dB)')
    ylabel('BER')
    title('64QAM环境下不同信道估计算法性能对比','Interpreter','none')

figure
    for i= 7:12
    eval(['semilogy(SNR_dB,ber_count_OTFS(i,:),''-',(dot(mod(i,6)+1,:)),''',''LineWidth'',2,''MarkerSize'',8)' ])
%      legend_str2{i}=['StudyCase ',Study_case,' with ',doppler_mode,' doppler in ',num2str(car_frequency/10^9),'GHz'];
    hold on
    end
%     for i= 1:3
%     eval(['semilogy(SNR_dB,ber_count_OFDM(i,:),''--',color(mod(i,4)+1,:),(dot(mod(i,4)+1,:)),''',''LineWidth'',2,''MarkerSize'',8)' ])
% %      legend_str2{i}=['StudyCase ',Study_case,' with ',doppler_mode,' doppler in ',num2str(car_frequency/10^9),'GHz'];
%     hold on
%     end
%     legend(legend_str2(9:16))
    axis([SNR_dB(1) SNR_dB(i_SNR) 1e-4 1e-0])
    grid on
    xlabel('SNR(dB)')
    ylabel('BER')
    title('64QAM环境下不同信道估计算法性能对比','Interpreter','none')

