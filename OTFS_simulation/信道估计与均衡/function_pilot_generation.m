function pilot_symbols = function_pilot_generation(QAM_mode,M_data,N_data,pilot_form,embedded_pilot_amp)
%FUNCTION_PILOT_GENERATION 此处显示有关此函数的摘要
%   此处显示详细说明

if strcmp(pilot_form, 'block')
    pilot_symbols=qammod(randi([0 QAM_mode-1],M_data,1),QAM_mode); % 生成导频
elseif strcmp(pilot_form, 'comb_like')
    pilot_symbols=qammod(randi([0 QAM_mode-1],1,N_data),QAM_mode); 
elseif strcmp(pilot_form, 'lattice')
    msg = '暂不支持格状导频';
    error(msg)
elseif strcmp(pilot_form, 'OTFS_delta')  
    pilot_symbols = zeros(1,N_data); 
    % pilot_symbols(ceil(N_data/2)) =qammod(randi([0 QAM_mode-1],1,1),QAM_mode);
    pilot_symbols(ceil(N_data/2)) =qammod(QAM_mode-1,QAM_mode,'bin')*embedded_pilot_amp;
else 
    msg = '请选择正确的导频形式：块状导频block，梳状导频comb_like，格状导频lattice，适用于嵌入式导频的OTFS_delta';
    error(msg)
    
end




end

