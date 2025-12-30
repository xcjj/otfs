clc
clear


% 生成 频域导频
M = 64;
N = 64;


pilot  = zeros(1,M) ;
pilot(1,M/2) = M;

pilot_delay = ifft(pilot) ;
plot(abs(pilot_delay));