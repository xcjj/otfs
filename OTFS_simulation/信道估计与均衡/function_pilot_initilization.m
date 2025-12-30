function [M_pilot,N_pilot,M_data,N_data,pilot_Indx,data_Indx] = function_pilot_initilization(M,N,channel_estimation_mode,pilot_form)
%FUNCTION_PILOT_INITILIZATION 此处显示有关此函数的摘要
%   此处显示详细说明

if (strcmp(channel_estimation_mode,'LS')||strcmp(channel_estimation_mode,'LMMSE')||strcmp(channel_estimation_mode,'ideal'))
% 9. 导频格式(块状导频block，格状导频lattice，梳状导频comb_like)
    if strcmp(pilot_form,'block')
        pilot_Inter=8; % 导频间隔 pilot interval,导频间隔是8意味着导频间计数需要+9
        % 计算导频和数据数量
        M_pilot = 0;
        N_pilot=ceil(N/(pilot_Inter+1)); % 向上取整
        if rem(N,pilot_Inter+1)==0 % 求余运算，如果被整除，则在末尾加一个导频
            N_pilot=N_pilot+1;
        end
        % 数据码元的数量
        M_data = M;
        N_data = N-N_pilot;
        % 导频与信息在一帧中位置的计算
        pilot_Indx=zeros(1,N_pilot);
        for i=1:N_pilot-1
            pilot_Indx(1,i)=(i-1)*(pilot_Inter+1)+1;
        end
        pilot_Indx(1,N_pilot)=N; % 默认最后一位必须是导频
        
        data_Indx = 1:N;
        data_Indx(pilot_Indx) =[];
    elseif strcmp(pilot_form,'comb_like')
        pilot_Inter=8; % 导频间隔 pilot interval,导频间隔是8意味着导频间计数需要+9
        % 计算导频和数据数量
        M_pilot=ceil(M/(pilot_Inter+1)); % 向下取整
        if rem(M,pilot_Inter+1)==0 % 求余运算，如果被整除，则在末尾加一个导频
            M_pilot=M_pilot+1;
        end
        N_pilot = 0;
        % 数据码元的数量
        N_data = N;
        M_data = M-M_pilot;
        % 导频与信息在一帧中位置的计算
        pilot_Indx=zeros(M_pilot,1);
        for i=1:M_pilot-1
            pilot_Indx(i,1)=(i-1)*(pilot_Inter+1)+1;
        end
        pilot_Indx(M_pilot,1)=M; % 默认最后一位必须是导频
        data_Indx = (1:M)';
        data_Indx(pilot_Indx) =[];
    elseif strcmp(pilot_form,'lattice')
        %%%%%%
        msg = '格状导频没做好。选择另外两个';
        error(msg)
    else
        msg = '导频格式参数输入错误，请重新检查参数,(如果是理想信道估计请把导频格式改成block)';
        error(msg) 
    end
elseif  strcmp(channel_estimation_mode,'embedded_pilot')
    if strcmp(pilot_form,'OTFS_delta')
    % 为表公平性，OTFS的DD域插入导频数量与其它信道估计插入导频数量相同
    % 计算方式沿用梳状导频
        pilot_Inter=8; % 导频间隔 pilot interval,导频间隔是8意味着导频间计数需要+9
        % 计算导频和数据数量
        M_pilot=ceil(M/(pilot_Inter+1)); % 向下取整
        if rem(M,pilot_Inter+1)==0 % 求余运算，如果被整除，则在末尾加一个导频
            M_pilot=M_pilot+1;
        end
        N_pilot = 0;
        % 数据码元的数量
        N_data = N;
        M_data = M-M_pilot;
        % 导频与信息在一帧中位置的计算
        
        pilot_Indx=(floor((M-M_pilot)/2+1):floor((M+M_pilot)/2))';
        data_Indx = (1:M)';
        data_Indx(pilot_Indx) =[];

    else
        msg = 'embedded_pilot所使用的导频格式与规定不符，不是OTFS_delta，请重新检查参数';
        error(msg)
    end
else 
    msg = '信道估计算法的变量输入错误，请重新检查参数';
    error(msg)
end
end

