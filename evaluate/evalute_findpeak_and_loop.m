clear; clc; close all;

%% === 1. 参数设置 ===
Fs = 200e6;              % 采样率 200MHz
T = 1e-4;                % 信号总时长 1ms
N_samples = round(Fs*T);
true_delay = 5;          % 真实时延 (采样点数)
win_len = 1024;          % 窗口长度
step_size = 256;         % 遍历步长 (窗口长度的1/4)
N_trials = 300;          % 重复实验次数
SNR = 10;                % 信噪比 (dB)

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
    sig_clean = zeros(1, N_samples);
    num_pulses = 50; 
    pulse_locs = sort(randi([win_len, N_samples-win_len], 1, num_pulses));
    
    for p = 1:num_pulses
        width = randi([5, 20]);
        p_start = max(1, pulse_locs(p)-width);
        p_end = min(N_samples, pulse_locs(p)+width);
        L = p_end - p_start + 1;
        sig_clean(p_start:p_end) = sig_clean(p_start:p_end) + (0.5 + 0.5*rand()) * hann(L)'; 
    end
    
    sig1 = sig_clean;
    sig2_clean = [zeros(1, true_delay), sig_clean(1:end-true_delay)];
    
    sig1 = awgn(sig1, SNR, 'measured');
    sig2 = awgn(sig2_clean, SNR, 'measured');
    
    % --- B. 方法1: 遍历窗口法 (Sliding Window) ---
    num_wins = floor((N_samples - win_len) / step_size);
    delays_sliding_raw = zeros(num_wins, 1);
    
    for w = 1:num_wins
        start_idx = (w-1)*step_size + 1;
        end_idx = start_idx + win_len - 1;
        
        chunk1 = sig1(start_idx:end_idx);
        chunk2 = sig2(start_idx:end_idx);
        
        [xc, lags] = xcorr(chunk1, chunk2);
        [~, max_idx] = max(abs(xc));
        delays_sliding_raw(w) = -lags(max_idx);
    end
    
    % === 【核心修改】 分组择优逻辑 ===
    % 物理含义：假设每 4 个连续窗口(覆盖1个win_len)中，必然有一个是对得比较准的
    % 我们取这4个里误差最小的那个作为这一组的代表
    group_size = 4; 
    num_groups = floor(num_wins / group_size);
    valid_len = num_groups * group_size;
    
    % 1. 计算所有窗口的原始误差
    raw_errors = abs(delays_sliding_raw(1:valid_len) - true_delay);
    
    % 2. 变形为矩阵 [4 x num_groups]
    errors_reshaped = reshape(raw_errors, group_size, num_groups);
    
    % 3. 列求最小 (择优)
    best_errors = min(errors_reshaped, [], 1);
    
    % 4. 记录本次实验的平均误差
    err_sliding_all(k) = mean(best_errors);
    
    
    % --- C. 方法2: 寻峰法 (Peak Finding) ---
    threshold = mean(sig1) + 3*std(sig1); 
    [pks, locs] = findpeaks(sig1, 'MinPeakHeight', threshold, 'MinPeakDistance', win_len/4);
    
    delays_peak = [];
    for p = 1:length(locs)
        center_loc = locs(p);
        if center_loc - win_len/2 < 1 || center_loc + win_len/2 > N_samples, continue; end
        idx_range = (center_loc - win_len/2) : (center_loc + win_len/2 - 1);
        
        [xc, lags] = xcorr(sig1(idx_range), sig2(idx_range));
        [~, max_idx] = max(abs(xc));
        delays_peak = [delays_peak; -lags(max_idx)];
    end
    
    if isempty(delays_peak)
        err_peak_all(k) = NaN; 
    else
        err_peak_all(k) = mean(abs(delays_peak - true_delay));
    end
    
    % 缓存最后一次数据用于绘图
    if k == N_trials
        last_sig1 = sig1;
        last_threshold = threshold;
        last_peak_locs = locs;
    end
end

% 数据清洗
valid_idx = ~isnan(err_peak_all);
err_peak_all = err_peak_all(valid_idx);
err_sliding_all = err_sliding_all(valid_idx);
N_valid = length(err_peak_all);

%% === 3. 绘图展示 (2x2 布局) ===
figure('Color', 'w', 'Position', [100, 100, 1200, 800]);

% --- 子图1 (左上): 遍历窗口示意图 (保持原样) ---
subplot(2, 2, 1);
plot(last_sig1, 'Color', [0.7 0.7 0.7]); hold on;
title('遍历窗口法示意', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('幅值'); xlabel('采样点');
xlim([0, 4096]); ylim([-1.5, 2]); 
for i = 1:4:16
    idx_s = (i-1)*step_size + 1;
    rectangle('Position', [idx_s, -1.2, win_len, 2.8], 'EdgeColor', 'b', 'LineStyle', '--', 'LineWidth', 1);
    text(idx_s, -1.4, 'Win', 'Color', 'b', 'FontSize', 8);
end
legend('含噪信号', '固定步长窗口', 'Location', 'northeast');
grid on;

% --- 子图2 (右上): 寻峰截取示意图 (保持原样) ---
subplot(2, 2, 2);
plot(last_sig1, 'Color', [0.7 0.7 0.7]); hold on;
title('寻峰截取法示意', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('幅值'); xlabel('采样点');
xlim([0, 4096]); ylim([-1.5, 2]);
yline(last_threshold, 'g--', 'Threshold');
valid_peaks_plot = last_peak_locs(last_peak_locs < 4096);
plot(valid_peaks_plot, last_sig1(valid_peaks_plot), 'rx', 'MarkerSize', 5, 'LineWidth', 2);
for i = 1:length(valid_peaks_plot)
    c_loc = valid_peaks_plot(i);
    w_start = c_loc - win_len/2;
    rectangle('Position', [w_start, -1.2, win_len, 2.8], 'EdgeColor', 'r', 'LineWidth', 1.5);
end
legend('含噪信号', '阈值线', '识别到的峰值', '自适应窗口', 'Location', 'northeast');
grid on;

% --- 子图3 (左下): 遍历法误差统计 (修正为分组择优后) ---
subplot(2, 2, 3);
plot(1:N_valid, err_sliding_all, '.-', 'Color', [0.2 0.4 0.8], 'LineWidth', 0.5); hold on;
mean_slide = mean(err_sliding_all);
yline(mean_slide, 'k--', 'LineWidth', 1.5);
title({'遍历法 - 300次实验误差', ['平均误差: ' num2str(mean_slide, '%.4f') ' (Samples)']}, ...
      'FontSize', 12, 'FontWeight', 'bold');
xlabel('实验序号'); ylabel('误差 (Samples)');
grid on;
ylim([0, 2]); % 锁定范围以便观察精度

% --- 子图4 (右下): 寻峰法误差统计 (保持原样) ---
subplot(2, 2, 4);
plot(1:N_valid, err_peak_all, '.-', 'Color', [0.8 0.2 0.2], 'LineWidth', 0.5); hold on;
mean_peak = mean(err_peak_all);
yline(mean_peak, 'k--', 'LineWidth', 1.5);
title({'寻峰法 - 300次实验误差', ['平均误差: ' num2str(mean_peak, '%.4f') ' (Samples)']}, ...
      'FontSize', 12, 'FontWeight', 'bold');
xlabel('实验序号'); ylabel('误差 (Samples)');
grid on;
ylim([0, 2]); % 锁定范围以便对比

% 主标题
sgtitle(['时延估计对比：遍历 vs 寻峰 (SNR=' num2str(SNR) 'dB)'], 'FontSize', 16);

fprintf('实验完成。\n');
fprintf('遍历法 平均误差: %.4f 点\n', mean_slide);
fprintf('寻峰法   平均误差: %.4f 点\n', mean_peak);