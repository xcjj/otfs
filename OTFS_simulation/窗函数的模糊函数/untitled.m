clc
clear

delta_f = 15*10^3;
delta_t = 1/delta_f;

L = 32;
M = 4;

sample_rate_t = 1; %采样倍率，10倍，使用时乘N
sample_rate_f = 1; %采样倍率，10倍，使用时乘M

t = 0:delta_t/(sample_rate_t*L):delta_t; % 时间轴
f = -delta_f*M:delta_f/(sample_rate_f):delta_f*M; % 频率轴


% 生成窗函数：

    % 单位能量的矩形窗
    Gtx =  1/sqrt(delta_t)*ones(size(t));
    Grx =  1/sqrt(delta_t)*ones(size(t));

expf = exp(-1i*2*pi*f'*t);
for i = 1:size(f,2)
[r(i,:),lags] = xcorr(Gtx,Grx.*expf(i,:));
end
[tt,ff] = meshgrid(lags/(sample_rate_t*L),f/(delta_f));





figure
surf(tt,ff,abs(r));
xlabel('时间t，-T到T')
ylabel('频率f，，-MΔf到MΔf，每一格代表子载波间隔Δf')

figure
for j = 1:sample_rate_f:sample_rate_f*M+1
 plot(lags/(sample_rate_t*L),real(r(j,:)));

JJ{round(j/sample_rate_f)+1} = ['f=' num2str(round(j/sample_rate_f)+1) 'Δf'];

hold on
end
legend(JJ)
xlabel('时间t，每一格代表周期T')

figure
for k = 1:sample_rate_t*8:sample_rate_t*L+1
% for k = sample_rate_t*L+1
plot(f/(delta_f),real(r(:,k)));

KK{round(k/(sample_rate_t*8))+1} = ['t=' num2str((k-1)/sample_rate_t/L) 'T'];

hold on
end
legend(KK)
xlabel('频率f，，-MΔf到MΔf，每一格代表子载波间隔Δf')

% r0 = r(:,sample_rate_t*L+1);
% % figure
% % % plot(f/(delta_f),r0);
% rtao = ifft(r0);
% figure
% plot(-M*sample_rate_f:M*sample_rate_f,abs(rtao));
% 
% rf0 = r0(sample_rate_f*M+10:sample_rate_f:sample_rate_f*2*M+1);
% figure
% plot(0:M-1,abs(rf0));
% 
% rtao0 = ifft(rf0);
% figure
% plot(0:M-1,abs(rtao0));
