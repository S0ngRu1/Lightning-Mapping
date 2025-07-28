%% ========================================================================
%  二维定位算法
% =========================================================================
clear; clc; close all;

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
pulses_per_group = 50; % 定义每组包含n个脉冲
% %从化局阈值
% noise = read_signal('..\\2024 822 85933.651462CH1.dat',1e8,1e8);
% filtered_noise = filter_bp(noise,30e6,80e6,5);
% threshold = mean(filtered_noise)+5*std(filtered_noise);
% --- 文件写入准备 ---
filename = 'result_yld_PULSE_CENTRIC_FINAL_'  + string(pulses_per_group) + '.txt';
fileID = fopen(filename, 'w');
fprintf(fileID, '%-13s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s\n', ...
    'Start_loc','Pulse_Len','t12', 't13', 't23', 'cos_alpha_opt', 'cos_beta_opt','Azimuth', 'Elevation', 'Rcorr', 't123');
num_total_blocks = numel(all_start_signal_loc)-1;

h_overall = waitbar(0, '正在初始化处理...', 'Name', '整体处理进度');
for j = 1:num_total_blocks

    progress = j / num_total_blocks;
    waitbar(progress, h_overall, sprintf('整体进度: 正在处理信号块 %d / %d', j, num_total_blocks));
    current_block_start = all_start_signal_loc(j);
    % --- 1. 读取当前处理块的完整信号 ---
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

    if isempty(pulse_catalog_in_block)
        fprintf('      在本块内未找到有效脉冲，跳过。\n');
        continue;
    end

    num_pulses_in_block = numel(pulse_catalog_in_block);
    fprintf('      在本块内找到 %d 个精确脉冲，开始处理...\n', num_pulses_in_block);

    % --- 遍历当前块内找到的【脉冲目录】 ---
    for pi = 1:pulses_per_group: num_pulses_in_block
        start_group_idx = pi;
        end_group_idx = min(pi + pulses_per_group - 1, num_pulses_in_block); % 防止末尾不足一组时越界
        pulse_group = pulse_catalog_in_block(start_group_idx:end_group_idx);
        % 如果组内脉冲太少，可以跳过
        if numel(pulse_group) < 2
            continue;
        end


        % **  确定该组的总起点和总终点 **
        group_start_idx = pulse_group(1).start_idx;
        group_end_idx   = pulse_group(end).end_idx;

        % --- 3. 根据总起止点，截取包含整组脉冲的长信号片段 ---
        signal1 = filtered_signal1(group_start_idx : group_end_idx);
        signal2 = filtered_signal2(group_start_idx : group_end_idx);
        signal3 = filtered_signal3(group_start_idx : group_end_idx);

        pulse_len = length(signal1); % 这是整个脉冲组的总长度

        % ** 长度对齐 (零填充)，以保证互相关等操作正常进行 **
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
        % 使用精确的亚采样点时刻作为定位点的时间戳
        ts_ns = 1/fs*1e9;
        absolute_sample_location = current_block_start + group_start_idx - 1;
        fprintf(fileID, '%-13d%-15d%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f\n', ...
            absolute_sample_location, pulse_len, t12, t13, t23, cos_alpha_opt, cos_beta_opt, Az_deg, El_deg, Rcorr,t123);
    end

end
close(h_overall);
fclose(fileID);