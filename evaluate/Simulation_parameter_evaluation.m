%% 雷电定位算法评估：独立绘图版 (上采样对比 & 自适应窗口)
% 采样率：200 MHz
clear; clc; close all;

%% === 1. 全局仿真参数设置 ===
fs = 200e6;                 
T_sim = 50e-6;              % 仿真时长 50us
N_total = round(T_sim * fs); 
true_delay_samples = 2.4;   % 亚采样延迟真值
num_trials = 300;           % 蒙特卡洛次数 
SNR_dB_range = -10:5:20;    % 信噪比范围

% 滤波器设计 (20M - 80M)
f1 = 20e6; f2 = 80e6;
[b_filt, a_filt] = butter(4, [f1, f2]/(fs/2)); 

fprintf('仿真开始 (Parfor并行加速中)...\n');

%% ========================================================================
%  实验一：上采样策略深度对比 (生成 Figure 1)
%  目的：证明 "原信号上采样" 优于 "互相关结果插值"
% ========================================================================
fprintf('正在进行实验一：上采样策略对比...\n');

% 固定参数：使用中等长度信号和匹配窗口
fixed_win = 1024;
fixed_sig_len = 800; 

% 待评估的上采样倍数
test_upsamples = [1, 5, 10, 30, 50]; 
corr_interp_compare = 8; % 对照组

% 1. 策略A：原信号上采样 
res_sig_up = zeros(length(test_upsamples), length(SNR_dB_range));
for i = 1:length(test_upsamples)
    curr_up = test_upsamples(i);
    fprintf('  -> [策略A] 原信号上采样 K=%d\n', curr_up);
    res_sig_up(i, :) = run_core_sim(SNR_dB_range, num_trials, fixed_win, ...
                                    curr_up, 0, ... % 0=无互相关插值
                                    fixed_sig_len, fs, N_total, true_delay_samples, b_filt, a_filt);
end

% 2. 策略B：互相关插值 (对照组)
fprintf('  -> [策略B] 互相关插值 M=%d (对照组)\n', corr_interp_compare);
res_corr_up = run_core_sim(SNR_dB_range, num_trials, fixed_win, ...
                           1, corr_interp_compare, ... % K=1, M=10
                           fixed_sig_len, fs, N_total, true_delay_samples, b_filt, a_filt);

%% ========================================================================
%  实验二：自适应窗口长度评估 (生成 Figure 2)
%  目的：证明不同密度的信号需要匹配不同的窗口
% ========================================================================
fprintf('正在进行实验二：自适应窗口匹配评估...\n');

% 定义4种信号场景 (模拟不同脉冲密度/持续时间)
scenarios = [400, 900, 1900, 3900];
scenario_names = {'高密度信号', '中等密度信号', ...
                  '低密度信号', '稀疏信号'};
% 待评估的窗口长度
test_wins = [512, 1024, 2048, 4096];
fixed_upsample = 10; % 固定一个合理的上采样率

res_exp2 = cell(1, 4); 

for s = 1:length(scenarios)
    curr_len = scenarios(s);
    fprintf('  -> 正在仿真场景 %d: 信号长度 %d\n', s, curr_len);
    
    scene_res = zeros(length(test_wins), length(SNR_dB_range));
    for w = 1:length(test_wins)
        curr_win = test_wins(w);
        scene_res(w, :) = run_core_sim(SNR_dB_range, num_trials, curr_win, ...
                                       fixed_upsample, 0, ...
                                       curr_len, fs, N_total, true_delay_samples, b_filt, a_filt);
    end
    res_exp2{s} = scene_res;
end

%% ========================================================================
%  绘图部分：生成两张独立图表
% ========================================================================

% --- Figure 1: 上采样策略对比 ---
figure('Units', 'pixels', 'Position', [100, 100, 800, 600], 'Color', 'w');
hold on; grid on; box on;

% 绘制策略A (彩色实线)
colors = jet(length(test_upsamples));
legend_str = {};
for i = 1:length(test_upsamples)
    plot(SNR_dB_range, res_sig_up(i, :), 'o-', 'LineWidth', 1.5, ...
         'Color', colors(i,:), 'MarkerFaceColor', colors(i,:), 'MarkerSize', 6);
    legend_str{end+1} = sprintf('原信号上采样 (K=%d)', test_upsamples(i));
end

% 绘制策略B (黑色虚线，加粗强调)
plot(SNR_dB_range, res_corr_up, 'k--s', 'LineWidth', 1.5, 'MarkerSize', 5, ...
     'MarkerFaceColor', 'k');
legend_str{end+1} = sprintf('对照组: 互相关插值 (M=%d)', corr_interp_compare);

title('不同上采样策略的时延估计精度对比', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('信噪比 SNR (dB)', 'FontSize', 12);
ylabel('时延估计标准差 (秒)', 'FontSize', 12);
legend(legend_str, 'Location', 'northeast', 'FontSize', 10);
set(gca, 'YScale', 'log', 'FontSize', 11); % 对数坐标展示数量级差异
ylim([1e-10, 1e-7]); % 根据需要调整范围


% --- Figure 2: 自适应窗口优势 ---
figure('Units', 'pixels', 'Position', [150, 150, 1000, 700], 'Color', 'w');
t = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

win_colors = lines(4); % 为4种窗口分配固定颜色

for s = 1:4
    nexttile;
    hold on; grid on; box on;
    data = res_exp2{s};
    
    for w = 1:4
        % 样式逻辑：最佳窗口实线加粗，其他窗口虚线变细
        lw = 1.5; l_style = '--'; marker = 'none'; alpha_val = 0.5;
        p = plot(SNR_dB_range, data(w, :), 'LineStyle', l_style, 'Marker', marker, ...
            'LineWidth', lw, 'Color', win_colors(w,:), 'MarkerSize', 4);
        p.Color(4) = alpha_val; % 设置透明度
    end
    
    title(scenario_names{s}, 'FontWeight', 'bold');
    xlabel('SNR (dB)'); 
    ylabel('STD (s)');
    set(gca, 'YScale', 'log');
    ylim([1e-10, 1e-6]);
    legend({'Win=512', 'Win=1024', 'Win=2048', 'Win=4096'}, ...
               'Location', 'southwest', 'FontSize', 9); 
end

title(t, '不同信号密度下的窗口长度性能评估', 'FontSize', 14, 'FontWeight', 'bold');

fprintf('所有绘图完成。\n');


%% ========================================================================
%  核心仿真函数 
% ========================================================================
function results = run_core_sim(SNRs, trials, win_len, sig_upsample, corr_upsample, sig_len, fs, N, true_delay, b, a)
    results = zeros(1, length(SNRs));
    
    % 预计算频域相移因子
    f_axis = [0:ceil(N/2)-1, -floor(N/2):-1] * (fs/N);
    phase_shift = exp(-1j * 2 * pi * f_axis * (true_delay/fs));
    
    for i = 1:length(SNRs)
        snr_val = SNRs(i);
        errors = zeros(1, trials);
        
        parfor k = 1:trials
            % 1. 信号生成
            raw_sig = randn(1, N);
            mask = zeros(1, N);
            center = floor(N/2);
            start_p = max(1, center - floor(sig_len/2));
            end_p = min(N, start_p + sig_len - 1);
            mask(start_p:end_p) = 1;
            sig = raw_sig .* mask;
            
            % 2. 亚采样延迟
            SIG = fft(sig);
            sig_delayed = ifft(SIG .* phase_shift, 'symmetric');
            
            % 3. 精确SNR缩放
            current_sig_power = sum(sig.^2) / sig_len; 
            if current_sig_power == 0, current_sig_power = 1; end
            target_sig_power = 1 * 10^(snr_val/10); 
            scale = sqrt(target_sig_power / current_sig_power);
            
            r1 = sig * scale + randn(1, N);
            r2 = sig_delayed * scale + randn(1, N);
            
            % 4. 滤波与截窗
            r1 = filter(b, a, r1);
            r2 = filter(b, a, r2);
            w_start = floor(N/2) - floor(win_len/2);
            idx = w_start : (w_start + win_len - 1);
            idx = idx(idx>0 & idx<=N);
            w1 = r1(idx);
            w2 = r2(idx);
            
            % 5. 策略分支
            if sig_upsample > 1
                x_old = 1:length(w1);
                x_new = linspace(1, length(w1), length(w1)*sig_upsample);
                w1_proc = interp1(x_old, w1, x_new, 'spline');
                w2_proc = interp1(x_old, w2, x_new, 'spline');
                lag_unit_scale = 1 / sig_upsample; 
            else
                w1_proc = w1; w2_proc = w2;
                lag_unit_scale = 1;
            end
            
            [cc, lags] = xcorr(w1_proc, w2_proc);
            
            if corr_upsample > 0
                cc_len = length(cc);
                c_old = 1:cc_len;
                c_new = linspace(1, cc_len, cc_len * corr_upsample);
                cc_final = interp1(c_old, cc, c_new, 'spline');
                lags_final = linspace(lags(1), lags(end), length(cc_final));
                lags_final = lags_final * lag_unit_scale;
            else
                cc_final = cc;
                lags_final = lags * lag_unit_scale;
            end
            
            % 6. 统一的时延估计 (三点抛物线)
            [~, max_idx] = max(cc_final);
            
            if max_idx > 1 && max_idx < length(cc_final)
                y1 = cc_final(max_idx-1); 
                y2 = cc_final(max_idx); 
                y3 = cc_final(max_idx+1);
                delta = 0.5 * (y1 - y3) / (y1 - 2*y2 + y3);
            else
                delta = 0;
            end
            
            grid_resolution = lags_final(2) - lags_final(1);
            est_delay_pts = lags_final(max_idx) + delta * grid_resolution;
            
            est_delay_sec = est_delay_pts / fs;
            errors(k) = abs(est_delay_sec - (true_delay/fs));
        end
        results(i) = std(errors);
    end
end