%% ========================================================================
%  二维定位算法 [四档自适应窗口最终版]
% =========================================================================
clear; clc; close all;

N = 3;
c = 0.299792458;
fs = 200e6;
step = 1e6;
upsampling_factor = 50;
start_signal_loc = 3.8e8;
end_signal_loc = 4.2e8;
signal_length = end_signal_loc - start_signal_loc;
all_start_signal_loc = start_signal_loc:step:end_signal_loc;
% % 从化局
% angle12 = -2.8381; angle13 = 50.3964; angle23 = 120.6568;
% d12 = 41.6496; d13 = 36.9015; d23 = 35.4481;
%引雷点
d12 = 24.9586; d13 = 34.9335; d23 = 24.9675;
angle12 = -110.8477; angle13 = -65.2405; angle23 = -19.6541;

%引雷点阈值
noise = read_signal('..\\20240822165932.6610CH1.dat',1e8,1e8);
filtered_noise = filter_bp(noise,30e6,80e6,5);
threshold = mean(filtered_noise)+5*std(filtered_noise);
% %从化局阈值
% noise = read_signal('..\\2024 822 85933.651462CH1.dat',1e8,1e8);
% filtered_noise = filter_bp(noise,30e6,80e6,5);
% threshold = mean(filtered_noise)+5*std(filtered_noise);
% 打开一个文本文件用于写入运行结果
fileID = fopen('result_yld_window_ADAPTIVE_1e6_factor4.txt', 'w');
fprintf(fileID, '%-13s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s\n', ...
    'Start_loc','Win_Len','t12', 't13', 't23', 'cos_alpha_opt', 'cos_beta_opt','Azimuth', 'Elevation', 'Rcorr', 't123');
for j = 1:numel(all_start_signal_loc)-1

    current_block_start = all_start_signal_loc(j);
    current_block_end = all_start_signal_loc(j+1);

    fprintf('>>>>>> 正在处理信号块: %d -- %d \n', current_block_start, current_block_end);

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
    scout_pulse_catalog = find_pulses_advanced(filtered_signal1, 3.545, fs, 4);
    pulse_count_in_chunk = numel(scout_pulse_catalog);
    % 调用决策函数，获得当前信号块应该使用的窗口长度
    dynamic_window_len = get_adaptive_window_length_4tier(pulse_count_in_chunk, step);
    fprintf('      本块密度: %d, 决策窗口: %d\n', pulse_count_in_chunk, dynamic_window_len);

    [peaks_in_block, locs_in_block] = findpeaks(filtered_signal1, 'MinPeakHeight', threshold, 'MinPeakDistance', dynamic_window_len/4);
    if isempty(locs_in_block)
        fprintf('      在本块内未找到有效脉冲，跳过。\n');
        continue;
    end

    num_peaks_in_block = numel(locs_in_block);
    h = waitbar(0, sprintf('正在处理块内 %d 个峰值...', num_peaks_in_block));

    for pi = 1:num_peaks_in_block
        waitbar(pi / num_peaks_in_block, h);
        idx = locs_in_block(pi);
        % 确保峰值不超出信号范围
        % 使用 dynamic_window_len 截取窗口信号
        win_start_idx = max(1, idx - floor(dynamic_window_len / 2) + 1);
        win_end_idx = min(step, idx + floor(dynamic_window_len / 2));

        % 截取窗口信号
        signal1 = filtered_signal1(win_start_idx:win_end_idx);
        signal2 = filtered_signal2(win_start_idx:win_end_idx);
        signal3 = filtered_signal3(win_start_idx:win_end_idx);
        % 去直流分量并应用窗函数
        [ch1_new, ch2_new, ch3_new] = deal(...
            real(windowsignal(detrend(signal1))), ...
            real(windowsignal(detrend(signal2))), ...
            real(windowsignal(detrend(signal3))));

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
        absolute_loc = current_block_start + idx;
        % 写入计算后的数据
        fprintf(fileID, '%-13d%-15d%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f\n', ...
            absolute_loc,dynamic_window_len, t12, t13, t23, cos_alpha_opt, cos_beta_opt, Az_deg, El_deg, Rcorr,t123);
    end
    close(h);
end
% 关闭文件
fclose(fileID);
