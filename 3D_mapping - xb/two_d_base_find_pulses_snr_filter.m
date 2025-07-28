%% ========================================================================
%  二维定位算法 [最终版：脉冲分组 + SNR筛选器]
% =========================================================================
clear; clc; close all;
% --- 0. 初始化和参数定义 ---

N = 3; c = 0.299792458; fs = 200e6; step = 1e4; upsampling_factor = 50;
start_signal_loc = 3.8e8; end_signal_loc = 4.2e8;
all_start_signal_loc = start_signal_loc:step:end_signal_loc;
% 从化局
% angle12 = -2.8381; angle13 = 50.3964; angle23 = 120.6568;
% d12 = 41.6496; d13 = 36.9015; d23 = 35.4481;
%引雷点
d12 = 24.9586; d13 = 34.9335; d23 = 24.9675;
angle12 = -110.8477; angle13 = -65.2405; angle23 = -19.6541;
%引雷点阈值
noise = read_signal('..\\20240822165932.6610CH1.dat',1e8,1e8);
filtered_noise = filter_bp(noise,30e6,80e6,5);
noise_std = std(filtered_noise);
threshold_factor = 3;      % find_pulses_advanced 的阈值因子
merge_gap_samples = 10;   % 脉冲融合的间隙阈值
pulses_per_group = 10; % 定义每组包含n个脉冲
% %从化局阈值
% noise = read_signal('..\\2024 822 85933.651462CH1.dat',1e8,1e8);
% filtered_noise = filter_bp(noise,30e6,80e6,5);
% threshold = mean(filtered_noise)+5*std(filtered_noise);

snr_threshold = 10; % SNR阈值：信噪比低于此值的脉冲组将被丢弃 
signal_band = [30e6, 80e6]; % 信号频带
noise_band = [0, 30e6;       % 第一段噪声区: 0 到 30 MHz
               80e6, 100e6];  % 第二段噪声区: 80 到 100 MHz
% *******************************************************************

% --- 文件写入准备 (增加一列用于记录SNR) ---
filename = 'result_yld_PULSE_GROUP_'  + string(pulses_per_group) + '_SNR_FILTERED_'+string(snr_threshold) +'.txt';
fileID = fopen(filename, 'w');
fprintf(fileID, '%-13s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s\n', ...
    'Start_loc','Pulse_Len','SNR','t12', 't13', 't23', 'cos_alpha_opt', 'cos_beta_opt','Azimuth', 'Elevation', 'Rcorr', 't123');

% --- 主循环开始---
num_total_blocks = numel(all_start_signal_loc)-1;
h_overall = waitbar(0, '正在初始化处理...', 'Name', '整体处理进度');

for j = 1:num_total_blocks
    progress = j / num_total_blocks;
    waitbar(progress, h_overall, sprintf('整体进度: 正在处理信号块 %d / %d', j, num_total_blocks));
    current_block_start = all_start_signal_loc(j);
    %     引雷点
    ch1 = read_signal('..\\20240822165932.6610CH1.dat', step, current_block_start);
    ch2 = read_signal('..\\20240822165932.6610CH2.dat', step, current_block_start);
    ch3 = read_signal('..\\20240822165932.6610CH3.dat', step, current_block_start);
    %     从化局
    % ch1 = read_signal('..\\2024 822 85933.651462CH1.dat',step,current_block_start);
    % ch2 = read_signal('..\\2024 822 85933.651462CH2.dat',step,current_block_start);
    % ch3 = read_signal('..\\2024 822 85933.651462CH3.dat',step,current_block_start+215/5);
    
    filtered_signal1 = filter_bp(ch1, 30e6, 80e6, 5);
    filtered_signal2 = filter_bp(ch2, 30e6, 80e6, 5);
    filtered_signal3 = filter_bp(ch3, 30e6, 80e6, 5);
    pulse_catalog_in_block = find_pulses_advanced(filtered_signal1, noise_std, fs, threshold_factor, merge_gap_samples);
    if isempty(pulse_catalog_in_block), continue; end
    num_pulses_in_block = numel(pulse_catalog_in_block);
    fprintf('      在本块内找到 %d 个精确脉冲，开始分组处理...\n', num_pulses_in_block);

    % --- 内循环：遍历脉冲组  ---
    for pi = 1:pulses_per_group:num_pulses_in_block
        start_group_idx = pi;
        end_group_idx = min(pi + pulses_per_group - 1, num_pulses_in_block);
        pulse_group = pulse_catalog_in_block(start_group_idx:end_group_idx);
        if numel(pulse_group) < 2, continue; end
        group_start_idx = pulse_group(1).start_idx;
        group_end_idx   = pulse_group(end).end_idx;
        signal1 = filtered_signal1(group_start_idx : group_end_idx);
        signal2 = filtered_signal2(group_start_idx : group_end_idx);
        signal3 = filtered_signal3(group_start_idx : group_end_idx);
        pulse_len = length(signal1);
        
        % 1. 计算当前脉冲组信号的带内信噪比 
        snr_value = calculate_in_band_snr(signal1, fs, signal_band, noise_band);
        
        % 2. 进行判断
        if snr_value < snr_threshold
            continue; % 跳过这个低质量的脉冲组
        end
        % --- 只有通过SNR检验的高质量信号，才能进入后续处理 ---
        max_len = max([length(signal1), length(signal2), length(signal3)]);
        signal1_padded = [signal1; zeros(max_len - length(signal1), 1)];
        signal2_padded = [signal2; zeros(max_len - length(signal2), 1)];
        signal3_padded = [signal3; zeros(max_len - length(signal3), 1)];

        % --- 后续的所有处理都使用填充对齐后的信号 ---
        [ch1_new, ch2_new, ch3_new] = deal(...
            real(windowsignal(detrend(signal1_padded))), ...
            real(windowsignal(detrend(signal2_padded))), ...
            real(windowsignal(detrend(signal3_padded))));


        % 上采样
        [ch1_up, ch2_up, ch3_up] = deal(...
            upsampling(ch1_new, upsampling_factor)', ...
            upsampling(ch2_new, upsampling_factor)', ...
            upsampling(ch3_new, upsampling_factor)');
        ch1_upsp = ch1_up(:,2);
        ch2_upsp = ch2_up(:,2);
        ch3_upsp = ch3_up(:,2);

        %互相关
        [r12_gcc,lags12_gcc] = xcorr(ch1_upsp ,ch2_upsp,'normalized');
        [r13_gcc,lags13_gcc] = xcorr(ch1_upsp ,ch3_upsp,'normalized');
        [r23_gcc,lags23_gcc] = xcorr(ch2_upsp ,ch3_upsp,'normalized');
        R12_gcc = max(r12_gcc);
        R13_gcc = max(r13_gcc);
        R23_gcc = max(r23_gcc);
        t12_gcc = cal_tau(r12_gcc,lags12_gcc');
        t13_gcc = cal_tau(r13_gcc,lags13_gcc');
        t23_gcc = cal_tau(r23_gcc,lags23_gcc');

        %
        %         %从化局
        %         t12 = t12_gcc *0.1;
        %         t13 = t13_gcc *0.1+1.600061;
        %         t23 = t23_gcc *0.1+1.600061;


        %引雷场
        t12 = t12_gcc *0.1;
        t13 = t13_gcc *0.1;
        t23 = t23_gcc *0.1;

        cos_beta_0 =((c*t13*d12*sind(angle12))-(c*t12*sind(angle13)*d13))/(d13*d12*sind(angle12-angle13)) ;
        cos_alpha_0 = ((c*t12)/d12-cos_beta_0*cosd(angle12))/sind(angle12);
        if abs(cos_beta_0)>1 || abs(cos_alpha_0)>1
            continue;
        end
        x0 = [cos_alpha_0,cos_beta_0];
        % 调用lsqnonlin函数进行优化
        options = optimoptions('lsqnonlin', 'MaxIter', 1000, 'TolFun', 1e-6);
        x = lsqnonlin(@(x) objective(x, t12, t13, t23,'yld'), x0, [-1 -1],[1 1], options);
        % 输出最优的cos(α)和cos(β)值
        cos_alpha_opt = x(1);
        cos_beta_opt = x(2);
        if abs(cos_alpha_opt)>1 || abs(cos_beta_opt)>1
            continue;
        end
        Az = atan2( cos_alpha_opt,cos_beta_opt);
        if abs(cos_beta_opt/cos(Az)) > 1
            continue;
        end
        El = acos( cos_beta_opt/cos(Az) );
        % 将弧度转换为角度
        Az_deg = rad2deg(Az);
        El_deg = rad2deg(El);
        if Az_deg < 0
            Az_deg = Az_deg + 360;
        end

        t123 = t12 + t23 - t13;
        Rcorr = (R12_gcc + R13_gcc + R23_gcc) / 3;
        % --- 写入结果 ---
        ts_ns = 1/fs*1e9;
        absolute_sample_location = current_block_start + group_start_idx - 1;
        fprintf(fileID, '%-13d%-15d%-15.2f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f\n', ...
            absolute_sample_location, pulse_len, snr_value, t12, t13, t23, cos_alpha_opt, cos_beta_opt, Az_deg, El_deg, Rcorr,t123);
    end
end
close(h_overall);
fclose(fileID);

%% --- 新版辅助函数：计算带内信噪比 (支持多段噪声区) ---
function snr = calculate_in_band_snr(signal_segment, fs, signal_band, noise_bands)
    % signal_segment: 输入的信号片段
    % fs: 采样率
    % signal_band: 信号频带, e.g., [30e6, 80e6]
    % noise_bands: 【新】噪声频带矩阵, e.g., [0, 30e6; 80e6, 100e6]
    
    nfft = 2^nextpow2(length(signal_segment));
    fft_result = fft(signal_segment, nfft);
    
    % 创建频率轴 (只取正频率部分)
    f_axis = fs*(0:(nfft/2))/nfft;
    
    % 找到信号频带对应的FFT索引
    signal_indices = find(f_axis >= signal_band(1) & f_axis <= signal_band(2));
    
    % *******************************************************************
    % ************************* 核心修改之处 ****************************
    % ******* 循环处理多段噪声区，并将索引合并 *******
    % *******************************************************************
    all_noise_indices = []; % 初始化一个空数组用于存储所有噪声索引
    for i = 1:size(noise_bands, 1)
        current_noise_band = noise_bands(i, :);
        noise_indices_segment = find(f_axis >= current_noise_band(1) & f_axis <= current_noise_band(2));
        all_noise_indices = [all_noise_indices, noise_indices_segment];
    end
    % *******************************************************************
    
    if isempty(signal_indices) || isempty(all_noise_indices)
        snr = -1; % 无法计算
        return;
    end
    
    % 计算信号和噪声频带内的平均功率
    power_signal = mean(abs(fft_result(signal_indices)).^2);
    power_noise = mean(abs(fft_result(all_noise_indices)).^2); % 在合并后的索引上计算
    
    if power_noise < 1e-12
        power_noise = 1e-12;
    end
    
    snr = power_signal / power_noise;
end
