function [info_tx_pilot_DT,info_tx_pilot_domain] = function_pilot_insert(M,N,Fn_M,Fn_N,Fn_M_pilot,Fn_N_pilot,info_tx_DT,pilot_symbols,pilot_domain,pilot_form,Data_Indx,pilot_Indx)
% 发送阶段的导频插入
%   此处显示详细说明
[M_data,N_data] = size(info_tx_DT);
M_pilot = M-M_data;
N_pilot = N-N_data;

if strcmp(pilot_domain, 'DT')
    % 不需要更改
    info_tx_domain = info_tx_DT;
elseif strcmp(pilot_domain, 'FT')
    info_tx_domain = Fn_M*info_tx_DT; % 信号转换为FT域
elseif strcmp(pilot_domain, 'DD')
    info_tx_domain = info_tx_DT*Fn_N; % 信号转换为DD域
elseif strcmp(pilot_domain, 'FD')
    info_tx_domain = Fn_M*info_tx_DT*Fn_N; % 信号转换为FD域   
else 
    msg = '请选择正确的导频插入域 （DT，FT, DD，FD）';
    error(msg)
end

info_tx_pilot_domain=zeros(M,N);

if strcmp(pilot_form, 'block')
    info_tx_pilot_domain(:,Data_Indx) = info_tx_domain;
    info_tx_pilot_domain(:,pilot_Indx) = repmat(pilot_symbols,1,N_pilot);
elseif strcmp(pilot_form, 'comb_like')
    info_tx_pilot_domain(Data_Indx,:) = info_tx_domain;
    info_tx_pilot_domain(pilot_Indx,:) = repmat(pilot_symbols,M_pilot,1);
elseif strcmp(pilot_form, 'lattice')
    msg = '暂不支持格状导频';
    error(msg)
elseif strcmp(pilot_form, 'OTFS_delta')
    info_tx_pilot_domain(Data_Indx,:) = info_tx_domain;
    info_tx_pilot_domain(pilot_Indx,:) = 0;
    info_tx_pilot_domain(pilot_Indx(ceil(M_pilot/2)),:) = pilot_symbols;
else 
    msg = '请选择正确的导频形式：块状导频block，梳状导频comb_like，格状导频lattice，或者OTFS_delta';
    error(msg)
    
end


if strcmp(pilot_domain, 'DT')
    % 不需要更改
    info_tx_pilot_DT = info_tx_pilot_domain;
elseif strcmp(pilot_domain, 'FT')
    info_tx_pilot_DT = Fn_M_pilot'*info_tx_pilot_domain; 
elseif strcmp(pilot_domain, 'DD')
    info_tx_pilot_DT = info_tx_pilot_domain*Fn_N_pilot';
elseif strcmp(pilot_domain, 'FD')
    info_tx_pilot_DT = Fn_M_pilot'*info_tx_pilot_domain*Fn_N_pilot'; 
end
            

end

