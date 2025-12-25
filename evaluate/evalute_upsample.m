clear; clc; close all;

%% === 1. 信号生成与参数设置 ===
Fs = 10e6;              % 采样率 10MHz
N = 64;                 % 信号长度
t = (0:N-1)/Fs;
true_delay_samples = 5.3; % 真实时延
interp_factor = 10;     % 插值倍数

% 生成原始脉冲信号 (高斯脉冲)
sig_ref = exp(-(t - 0.4*max(t)).^2 / (2 * (0.08*max(t))^2));

% 生成带时延的信号 (利用频域相移生成精确的小数时延)
X = fft(sig_ref);
f = [0:N/2, -N/2+1:-1] * (Fs/N);
phase_shift = exp(-1j * 2 * pi * f * (true_delay_samples / Fs));
sig_delayed = ifft(X .* phase_shift, 'symmetric');

% 添加噪声
SNR = 60;
sig_ref = awgn(sig_ref, SNR, 'measured');
sig_delayed = awgn(sig_delayed, SNR, 'measured');

%% === 2. 方法一：信号插值法 (Signal Interpolation) ===
% 1. 信号插值
sig_ref_up = interp(sig_ref, interp_factor);
sig_del_up = interp(sig_delayed, interp_factor);

% 2. 互相关 (【修改点】：添加 'normalized' 参数)
[xc_1, lags_1] = xcorr(sig_ref_up, sig_del_up, 'normalized');

% 3. 寻峰与拟合
[~, max_idx_1] = max(xc_1);
y_pts_1 = xc_1(max_idx_1-1 : max_idx_1+1);
x_pts_1 = lags_1(max_idx_1-1 : max_idx_1+1);
[p_1, est_peak_x_1, est_peak_y_1] = parabolic_fit(x_pts_1, y_pts_1);
delay_est_1 = -est_peak_x_1 / interp_factor; 

%% === 3. 方法二：互相关插值法 (Xcorr Interpolation) ===
% 1. 原始互相关 (【修改点】：添加 'normalized' 参数)
[xc_raw, lags_raw] = xcorr(sig_ref, sig_delayed, 'normalized');

% 2. 互相关函数插值
xc_2 = interp(xc_raw, interp_factor);
% 重新生成插值后的 lags 轴
lags_2 = linspace(min(lags_raw), max(lags_raw), length(xc_2));

% 3. 寻峰与拟合
[~, max_idx_2] = max(xc_2);
y_pts_2 = xc_2(max_idx_2-1 : max_idx_2+1);
x_pts_2 = lags_2(max_idx_2-1 : max_idx_2+1);
[p_2, est_peak_x_2, est_peak_y_2] = parabolic_fit(x_pts_2, y_pts_2);
delay_est_2 = -est_peak_x_2; 

%% === 4. 绘图展示 ===
figure('Color', 'w', 'Position', [100, 50, 1000, 800]);

% --- 子图1 (顶部)：原始信号 vs 上采样信号 ---
subplot(2, 2, [1, 2]); 
t_up = linspace(min(t), max(t), length(sig_ref_up));
% 画原始信号
stem(t*1e6, sig_ref, 'filled', 'MarkerSize', 5, 'Color', [0.2 0.2 0.2], 'DisplayName', '原始采样点 (Fs)');
hold on;
% 画上采样信号
plot(t_up*1e6, sig_ref_up, 'r-', 'LineWidth', 1.5, 'DisplayName', ['上采样插值 (' num2str(interp_factor) 'x Fs)']);

title('原始信号 vs 上采样插值信号', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Time (\mus)'); ylabel('Amplitude');
legend('Location', 'northeast'); grid on;
xlim([min(t)*1e6, max(t)*1e6]); 

% --- 子图2 (左下)：信号插值法结果 ---
subplot(2, 2, 3);
display_range_1 = max_idx_1-30 : max_idx_1+30;
plot(lags_1(display_range_1)/interp_factor, xc_1(display_range_1), 'b.-', 'MarkerSize', 8); hold on;
% 拟合曲线
x_fit_1 = linspace(min(x_pts_1), max(x_pts_1), 100);
plot(x_fit_1/interp_factor, polyval(p_1, x_fit_1), 'r-', 'LineWidth', 2);
% 标记
plot(est_peak_x_1/interp_factor, est_peak_y_1, 'rx', 'MarkerSize', 12, 'LineWidth', 2);
xline(-true_delay_samples, 'g--', 'Label', 'True');

title({['信号预插值  (真实时延 = ' num2str(true_delay_samples) ' 点)'], ['Est: ' num2str(delay_est_1, '%.4f')]}, 'FontSize', 11, 'FontWeight', 'bold');
xlabel('Lags (Samples)'); ylabel('Correlation (Normalized)'); grid on;
legend('插值后互相关', '抛物线拟合', '估计峰值', 'Location', 'southwest');

% --- 子图3 (右下)：互相关插值法结果 ---
subplot(2, 2, 4);
display_range_2 = max_idx_2-30 : max_idx_2+30;
plot(lags_2(display_range_2), xc_2(display_range_2), 'k.-', 'MarkerSize', 8); hold on;
% 拟合曲线
x_fit_2 = linspace(min(x_pts_2), max(x_pts_2), 100);
plot(x_fit_2, polyval(p_2, x_fit_2), 'r-', 'LineWidth', 2);
% 标记
plot(est_peak_x_2, est_peak_y_2, 'rx', 'MarkerSize', 12, 'LineWidth', 2);
xline(-true_delay_samples, 'g--', 'Label', 'True');

title({['互相关后插值  (真实时延 = ' num2str(true_delay_samples) ' 点)'], ['Est: ' num2str(delay_est_2, '%.4f')]}, 'FontSize', 11, 'FontWeight', 'bold');
xlabel('Lags (Samples)'); ylabel('Correlation '); grid on;
legend('互相关插值结果', '抛物线拟合', '估计峰值', 'Location', 'southwest');

sgtitle(['时延估计流程与精度对比 (真实时延 = ' num2str(true_delay_samples) ' 点)'], 'FontSize', 14);

%% === 辅助函数 ===
function [p, peak_x, peak_y] = parabolic_fit(x, y)
    if iscolumn(x), x = x'; end
    if iscolumn(y), y = y'; end
    p = polyfit(x, y, 2);
    peak_x = -p(2) / (2 * p(1));
    peak_y = polyval(p, peak_x);
end