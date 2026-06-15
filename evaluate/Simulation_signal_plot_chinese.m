%% JGR 出版物信号处理工作流可视化
% 更新日志：
% 1. 减少至 2 个子图。
% 2. 使用平铺布局的共享 Y 轴标签。
% 3. 全中文标注。

clear; clc; close all;

%% === 1. 信号生成 ===
fs = 200e6;                     
N = 4096;                       
win_len = 1024;                 
snr_val = 15;                   
true_delay_samples = 20.5;      
compress_scale = 0.4; % 原始信号可视化缩放因子

% 滤波器设计：20-80 MHz 带通
f_low = 20e6; f_high = 80e6;

% --- 生成双极性脉冲群 ---
rng(42); 
sig = zeros(1, N);
num_pulses = 40; 
range_start = floor(N/2) - 400;
range_end = floor(N/2) + 400;
pulse_locs = sort(randi([range_start, range_end], 1, num_pulses));

for p = 1:num_pulses
    width = randi([5, 12]);
    pulse_sign = sign(randn()); 
    amp = pulse_sign * (0.5 + 0.5*rand());
    p_start = max(1, pulse_locs(p)-width);
    p_end = min(N, pulse_locs(p)+width);
    L = p_end - p_start + 1;
    sig(p_start:p_end) = sig(p_start:p_end) + amp * hann(L)'; 
end
sig = sig / max(abs(sig));

% --- 应用延迟与噪声 ---
f_axis = [0:ceil(N/2)-1, -floor(N/2):-1] * (fs/N);
phase_shift = exp(-1j * 2 * pi * f_axis * (true_delay_samples/fs));
SIG = fft(sig);
sig_delayed = ifft(SIG .* phase_shift, 'symmetric');
sig_power = mean(sig.^2);
noise_power = sig_power / (10^(snr_val/10));
noise_std = sqrt(noise_power);

% 生成原始含噪信号 (原始与延迟)
r1 = sig + noise_std * randn(1, N);
r2 = sig_delayed + noise_std * randn(1, N);

% --- 滤波 (零相位) 与加窗 ---
r1_filt = filter_bp(r1, f_low, f_high, 5);
r2_filt = filter_bp(r2, f_low, f_high, 5);
center_idx = floor(N/2);
w_start = center_idx - floor(win_len/2);
idx = w_start : (w_start + win_len - 1);
ham_win = hamming(win_len)';
w1_final = r1_filt(idx) .* ham_win;

%% === 2. 绘图 (中文版) ===
% 调整图形比例
figure('Units', 'inches', 'Color', 'w'); 
t = tiledlayout(2, 1, 'TileSpacing', 'tight', 'Padding', 'compact');

% 定义显示范围
plot_range = (w_start - 200) : (w_start + win_len + 200);
plot_range = plot_range(plot_range > 0 & plot_range <= N);

% 中文字体设置 (建议使用微软雅黑)
fontName = 'Microsoft YaHei'; 
fontSize = 10;
lineWidthMain = 1.2;
lineWidthThin = 0.8;

% === 面板 (a): 原始信号 vs 延迟信号 (原始数据) ===
ax1 = nexttile; hold on; box on;
plot(plot_range, r1(plot_range) * compress_scale, 'Color', [0.7 0.7 0.7], 'LineWidth', lineWidthThin, ...
    'DisplayName', '原始信号');
plot(plot_range, r2(plot_range) * compress_scale, 'Color', [0.85 0.325 0.098], 'LineWidth', lineWidthThin, ...
    'DisplayName', '延迟信号');

xlim([plot_range(1), plot_range(end)]);
ylim([-0.5, 0.5]); 
set(gca, 'XTickLabel', []); 
grid on; set(gca, 'GridLineStyle', ':');
legend('Location', 'northeast', 'Box', 'off', 'FontName', fontName);
text(0.02, 0.9, '(a)', 'Units', 'normalized', 'FontName', fontName, 'FontSize', 12, 'FontWeight', 'bold');

% === 面板 (b): 滤波与加窗处理过程 ===
ax2 = nexttile; hold on; box on;
% 滤波后 (未加窗背景)
plot(plot_range, r1_filt(plot_range), '--', 'Color', [0.3010 0.7450 0.9330], 'LineWidth', 1, ...
    'DisplayName', '带通滤波后');
% 加窗后 (最终结果)
plot(idx, w1_final, 'Color', [0 0 0.8], 'LineWidth', 1.5, ...
    'DisplayName', '最终加窗信号');
% 海明窗可视化
plot(idx, ham_win * max(abs(w1_final))*1.2, 'Color', [0.4660 0.6740 0.1880], 'LineStyle', ':', 'LineWidth', 1.5, ...
    'DisplayName', 'hann窗函数');

xlabel('时间 (采样点)', 'FontName', fontName, 'FontSize', fontSize);
xlim([plot_range(1), plot_range(end)]);
ylim([-0.1, 0.2]); 
grid on; set(gca, 'GridLineStyle', ':');

% 窗口边界线
xline(idx(1), 'k-', 'LineWidth', 0.5, 'Color', [0.4 0.4 0.4], 'HandleVisibility', 'off');
xline(idx(end), 'k-', 'LineWidth', 0.5, 'Color', [0.4 0.4 0.4], 'HandleVisibility', 'off');

legend('Location', 'northeast', 'Box', 'off', 'FontName', fontName);
text(0.02, 0.9, '(b)', 'Units', 'normalized', 'FontName', fontName, 'FontSize', 12, 'FontWeight', 'bold');

% === 添加共享 Y 轴标题 ===
ylabel(t, '幅度', 'FontName', fontName, 'FontSize', fontSize, 'FontWeight', 'normal');

% 全局字体统一调整
set(findall(gcf,'-property','FontSize'),'FontSize', fontSize);
set(findall(gcf,'-property','FontName'),'FontName', fontName);

fprintf('图表更新完成：已切换为中文显示并添加共享 Y 轴标签。\n');

%% 辅助函数
function filtered_signal = filter_bp(signal,f1,f2,order)
    Fs = 200e6;
    fn = Fs/2;
    Wn = [f1 f2]/fn;
    [b,a] = butter(order,Wn);
    filtered_signal = filtfilt(b,a,signal);
end