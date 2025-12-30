% Channel_DT=function_Gen_delay_time_channel(M,N,delay_taps,doppler_taps,chan_amplitude_coef_taps,length_CP,interval_mode); 

function [Channel_DT]=function_Gen_delay_time_channel(M,N,delay_taps,doppler_taps,chan_amplitude_coef_taps,length_interval,interval_mode)

z=exp(1i*2*pi/N/M);
taps=length(delay_taps);
delay_max=max(delay_taps);
if strcmp(interval_mode,'CP')
    length_CP = length_interval;
%     frame_size=(length_CP+M)*N;
    frame_size=(length_CP+M)*(N);
elseif strcmp(interval_mode,'ZP')
    length_ZP = length_interval;
    frame_size=(M+length_ZP)*N;
elseif strcmp(interval_mode,'RCP')    
    length_RCP = length_interval;
    frame_size=length_RCP+M*N;
elseif strcmp(interval_mode,'RZP')    
    length_RZP = length_interval;
    frame_size=M*N+length_RZP;
else
    frame_size=M*N;
end

Channel_DT=zeros(delay_max+1,frame_size);
for q=1:frame_size
    for taps_i=1:taps
        gain_i=chan_amplitude_coef_taps(taps_i);
        delay_taps_i=delay_taps(taps_i);
        doppler_i=doppler_taps(taps_i);  
        Channel_DT(delay_taps_i+1,q)=Channel_DT(delay_taps_i+1,q)+gain_i*z^((doppler_i)*(q-1-delay_taps_i));  
        % Channel_DT(delay_taps_i+1,q)=Channel_DT(delay_taps_i+1,q)+gain_i*z^((doppler_i)*(ceil(q/ceil(frame_size/N))*ceil(frame_size/N)-1-delay_taps_i));  
    end    
end
end