%https://ww2.mathworks.cn/help/comm/ref/comm.rayleighchannel-system-object.html?s_tid=srchtitle_comm.RayleighChannel_6

clear
clc

pskModulator = comm.PSKModulator;
insig = randi([0,pskModulator.ModulationOrder-1],1024,1);
channelInput = pskModulator(insig);

channelInput  = zeros(1024,1);
channelInput(1,1) = 1;
channelInput(20,1) = 1;


rayleighchan = comm.RayleighChannel( ...
    'SampleRate',10e3, ...
    'PathDelays',[0 1.5e-4], ...
    'AveragePathGains',[2 3], ...
    'NormalizePathGains',true, ...
    'MaximumDopplerShift',30, ...
    'DopplerSpectrum',{doppler('Gaussian',0.6),doppler('Flat')}, ...
    'RandomStream','mt19937ar with seed', ...
    'Seed',22, ...
    'PathGainsOutputPort',true);

[chanOut1,pathGains1] = rayleighchan(channelInput);

abs(pathGains1);



% release(rayleighchan);
% rayleighchan.RandomStream = 'Global stream';
% 
% rng(22)
% 
% [chanOut2,pathGains2] = rayleighchan(channelInput);
% 
% isequal(chanOut1,chanOut2)
% isequal(pathGains1,pathGains2)