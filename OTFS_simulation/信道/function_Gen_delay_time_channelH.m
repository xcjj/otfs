% Channel_DT=function_Gen_delay_time_channel(M,N,delay_taps,doppler_taps,chan_amplitude_coef_taps,length_CP,interval_mode); 

function [ChannelH_DT]=function_Gen_delay_time_channelH(M,N,delay_taps,doppler_taps, chan_amplitude_coef_taps,length_interval,interval_mode)

% interval_mode = 'None';
% interval_mode = 'CP';

z=exp(1i*2*pi/N/M);
taps=length(delay_taps);
delay_max=max(delay_taps);
if strcmp(interval_mode,'CP')
    length_CP = length_interval;
    frame_size_M=(length_CP+M);
    frame_size_N=N;
elseif strcmp(interval_mode,'ZP')
    length_ZP = length_interval;
    frame_size_M=(M+length_ZP);
    frame_size_N=N;
elseif strcmp(interval_mode,'RCP')    
    msg = '不支持RCP，RZP';
    error(msg)
%     length_RCP = length_interval;
%    frame_size=length_RCP+M*N;
elseif strcmp(interval_mode,'RZP')    
%     length_RZP = length_interval;
%     frame_size=M*N+length_RZP;
    msg = '不支持RCP，RZP';
    error(msg)
else
    frame_size_M=M;
    frame_size_N=N;
end

ChannelH_DT=zeros(frame_size_M,frame_size_N);
% for q=1:frame_size
%     for taps_i=1:taps
%         gain_i=chan_amplitude_coef_taps(taps_i);
%         delay_taps_i=delay_taps(taps_i);
%         doppler_i=doppler_taps(taps_i);  
%         Channel_DT(delay_taps_i+1,q)=Channel_DT(delay_taps_i+1,q)+gain_i*z^(doppler_i*(q-1-delay_taps_i));  
%     end    
% end


for taps_i=1:taps
        gain_i=chan_amplitude_coef_taps(taps_i);
        delay_taps_i=delay_taps(taps_i);
        doppler_i=doppler_taps(taps_i); 
%         if strcmp(doppler_mode,'time varying')
%             doppler_rate_i = doppler_rate_taps(taps_i);
%         elseif strcmp(doppler_mode,'fixed')
%             doppler_rate_i = 0;
%         end
        doppler_rate_i = 0;
    for t_i=1:frame_size_N
        t = (t_i-1)*M+delay_taps_i+1;
        ChannelH_DT(delay_taps_i+1,t_i)=ChannelH_DT(delay_taps_i+1,t_i)+gain_i*z^((doppler_i+doppler_rate_i*t*1.04*10^-6)*(t-delay_taps_i+1));  
    end    
end



end