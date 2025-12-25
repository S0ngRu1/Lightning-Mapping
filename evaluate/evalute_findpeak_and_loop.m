clear; clc; close all;

%% === 1. 参数设置 ===
Fs = 200e6;              % 采样率 200MHz
T = 1e-4;               % 信号总时长 1ms
N_samples = round(Fs*T);
true_delay = 5;        % 真实时延 (采样点数)
win_len = 1024;          % 窗口长度
step_size = 256;        % 遍历步长
N_trials = 300;         % 重复实验次数
SNR = 10;               % 信噪比 (dB)

% 存储每次实验的误差
err_sliding_all = zeros(N_trials, 1);
err_peak_all = zeros(N_trials, 1);

% 用于最后画示意图的缓存变量
last_sig1 = [];
last_threshold = 0;
last_peak_locs = [];

%% === 2. 循环实验 (300次) ===
fprintf('正在进行300次蒙特卡洛模拟...\n');

for k = 1:N_trials
    % --- A. 生成信号 ---
    % 模拟稀疏的闪电脉冲信号
    sig_clean = zeros(1, N_samples);
    num_pulses = 50; 
    
    % 随机生成脉冲位置
    pulse_locs = sort(randi([win_len, N_samples-win_len], 1, num_pulses));
    
    for p = 1:num_pulses
        width = randi([5, 20]);
        % 边界保护
        p_start = max(1, pulse_locs(p)-width);
        p_end = min(N_samples, pulse_locs(p)+width);
        L = p_end - p_start + 1;
        % 叠加脉冲
        sig_clean(p_start:p_end) = sig_clean(p_start:p_end) + (0.5 + 0.5*rand()) * hann(L)'; 
    end
    
    % 信号1 (基准) 和 信号2 (时延 + 噪声)
    sig1 = sig_clean;
    sig2_clean = [zeros(1, true_delay), sig_clean(1:end-true_delay)];
    
    % 添加噪声
    sig1 = awgn(sig1, SNR, 'measured');
    sig2 = awgn(sig2_clean, SNR, 'measured');
    
    % --- B. 方法1: 遍历窗口法 (Sliding Window) ---
    delays_sliding = [];
    num_wins = floor((N_samples - win_len) / step_size);
    
    for w = 1:num_wins
        start_idx = (w-1)*step_size + 1;
        end_idx = start_idx + win_len - 1;
        
        chunk1 = sig1(start_idx:end_idx);
        chunk2 = sig2(start_idx:end_idx);
        
        [xc, lags] = xcorr(chunk1, chunk2);
        [~, max_idx] = max(abs(xc));
        estimated_delay = -lags(max_idx);
        
        delays_sliding = [delays_sliding; estimated_delay];
    end
    err_sliding_all(k) = mean(abs(delays_sliding - true_delay));
    
    % --- C. 方法2: 寻峰法 (Peak Finding) ---
    threshold = mean(sig1) + 3*std(sig1); % 阈值设高一点，模拟只找显著信号
    [pks, locs] = findpeaks(sig1, 'MinPeakHeight', threshold, 'MinPeakDistance', win_len/4);
    
    delays_peak = [];
    for p = 1:length(locs)
        center_loc = locs(p);
        if center_loc - win_len/2 < 1 || center_loc + win_len/2 > N_samples
            continue;
        end
        idx_range = (center_loc - win_len/2) : (center_loc + win_len/2 - 1);
        chunk1 = sig1(idx_range);
        chunk2 = sig2(idx_range);
        
        [xc, lags] = xcorr(chunk1, chunk2);
        [~, max_idx] = max(abs(xc));
        estimated_delay = -lags(max_idx);
        
        delays_peak = [delays_peak; estimated_delay];
    end
    
    if isempty(delays_peak)
        err_peak_all(k) = NaN; 
    else
        err_peak_all(k) = mean(abs(delays_peak - true_delay));
    end
    
    % 保存最后一次的数据用于画示意图
    if k == N_trials
        last_sig1 = sig1;
        last_threshold = threshold;
        last_peak_locs = locs;
    end
end

% 移除无效数据
valid_idx = ~isnan(err_peak_all);
err_peak_all = err_peak_all(valid_idx);
err_sliding_all = err_sliding_all(valid_idx);
N_valid = length(err_peak_all);

%% === 3. 绘图展示 (2x2 布局) ===
figure('Color', 'w', 'Position', [100, 100, 1200, 800]);

% --- 子图1 (左上): 遍历窗口示意图 ---
subplot(2, 2, 1);
plot(last_sig1, 'Color', [0.7 0.7 0.7]); hold on;
title('遍历窗口法示意', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('幅值'); xlabel('采样点');
xlim([0, 4096]); ylim([-1.5, 2]); % 只展示前4096点，方便看清

for i = 1:4:16
    idx_s = (i-1)*step_size + 1;
    % 画虚线框
    rectangle('Position', [idx_s, -1.2, win_len, 2.8], 'EdgeColor', 'b', 'LineStyle', '--', 'LineWidth', 1);
    text(idx_s, -1.4, 'Win', 'Color', 'b', 'FontSize', 8);
end
legend('含噪信号', '固定步长窗口', 'Location', 'northeast');
grid on;

% --- 子图2 (右上): 寻峰截取示意图 ---
subplot(2, 2, 2);
plot(last_sig1, 'Color', [0.7 0.7 0.7]); hold on;
title('寻峰截取法示意', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('幅值'); xlabel('采样点');
xlim([0, 4096]); ylim([-1.5, 2]);
% 画阈值线
yline(last_threshold, 'g--', 'Threshold');
% 标出峰值点
valid_peaks_plot = last_peak_locs(last_peak_locs < 4096);
plot(valid_peaks_plot, last_sig1(valid_peaks_plot), 'rx', 'MarkerSize', 5, 'LineWidth', 2);
% 画寻峰框
for i = 1:length(valid_peaks_plot)
    c_loc = valid_peaks_plot(i);
    w_start = c_loc - win_len/2;
    % 画实线框
    rectangle('Position', [w_start, -1.2, win_len, 2.8], 'EdgeColor', 'r', 'LineWidth', 1.5);
    text(w_start, 1.4, 'Win', 'Color', 'r', 'FontSize', 8);
end
legend('含噪信号', '阈值线', '识别到的峰值', '自适应窗口', 'Location', 'northeast');
grid on;

% --- 子图3 (左下): 遍历法误差统计 ---
subplot(2, 2, 3);
plot(1:N_valid, err_sliding_all, '.-', 'Color', [0.2 0.4 0.8], 'LineWidth', 0.5); hold on;
mean_slide = mean(err_sliding_all);
title({'遍历法 - 300次实验误差', ['平均误差: ' num2str(mean_slide, '%.2f') ' (Samples)']}, 'FontSize', 12, 'FontWeight', 'bold');xlabel('实验序号'); ylabel('平均时延误差 (Samples)');
grid on;
ylim([0, max(err_sliding_all)*1.2]); % 自动调整范围

% --- 子图4 (右下): 寻峰法误差统计 ---
subplot(2, 2, 4);
plot(1:N_valid, err_peak_all, '.-', 'Color', [0.8 0.2 0.2], 'LineWidth', 0.5); hold on;
mean_peak = mean(err_peak_all);
title({'寻峰法 - 300次实验误差', ['平均误差: ' num2str(mean_peak, '%.2f') ' (Samples)']}, 'FontSize', 12, 'FontWeight', 'bold');xlabel('实验序号'); ylabel('平均时延误差 ');
grid on;
% 统一如果误差很小，设置一个固定的小范围以便观察微小波动，或者跟随左图量级对比
% 这里为了看清微小误差，自动缩放，但在标题里能看出数量级差异
if max(err_peak_all) < 2
    ylim([0, 2]);
else
    ylim([0, max(err_peak_all)*1.2]);
end

% 主标题
sgtitle(['时延估计方法对比：遍历窗口 vs 寻峰截取 (SNR=' num2str(SNR) 'dB)'], 'FontSize', 16);

% 输出文字结论
fprintf('实验完成。\n');
fprintf('遍历法平均误差: %.4f 点\n', mean_slide);
fprintf('寻峰法平均误差: %.4f 点\n', mean_peak);