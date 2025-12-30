% 信号在通过信道之后的时域表示
% 公式见电子书[R1] P63 4.4.1 公式（4.36）
% [R1]. Y. Hong, T. Thaj, E. Viterbo, ``Delay-Doppler Communications: Principles and Applications'', Academic Press, 2022, ISBN:9780323850285

function [r,noise]=function_Gen_channel_output_DT(channel_delay_set,Channel_DT,info_tx_P2S_DT,sigma_2)
r=zeros(size(info_tx_P2S_DT));
noise= sqrt(sigma_2/2)*(randn(size(info_tx_P2S_DT)) + 1i*randn(size(info_tx_P2S_DT)));
for q=0:(size(r,1)-1)
    for l=channel_delay_set
        if(q>=l)
            r(q+1)=r(q+1)+Channel_DT(l+1,q+1)*info_tx_P2S_DT(q-l+1);    
        end
    end
end

%  r=r+noise;
end

