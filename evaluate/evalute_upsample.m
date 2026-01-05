%% Figure 3 Generation: Spline Upsampling Schematic for JGR
% Goal: Illustrate the concept of Section 3.2.2 (Spline Kernel Upsampling)
% Format: JGR Standard (Arial, English labels, (a)/(b) panels)

clear; clc; close all;

%% === 1. Parameter Setup ===
fs_coarse = 20e6;           % 原始采样率 20 MHz (模拟粗粒度)
K = 10;                     % 上采样倍数 (Upsampling factor)
fs_fine = fs_coarse * K;    % 高倍采样率
T_total = 1e-6;             % 信号持续时间
t_coarse = 0:1/fs_coarse:T_total;
t_fine = 0:1/fs_fine:T_total;

% 生成一个具有亚采样时延的模拟信号 (Sinc-like pulse)
true_delay = 0.35 * (1/fs_coarse); % 设置一个非整数采样点的时延 (0.35个采样周期)
center_t = T_total/2;

% 定义连续时间波形函数 (用于画"真值"和生成采样点)
sig_fun = @(t, delay) exp(-(t - center_t - delay).^2 / (0.1e-6)^2) .* ...
                      cos(2*pi*5e6*(t - center_t - delay));

% 1. 生成数据
sig_true_continuous = sig_fun(t_fine, true_delay);   % "真实"的模拟信号 (用于背景参考)
sig_discrete = sig_fun(t_coarse, true_delay);        % 原始离散采样信号

% 2. 核心算法：三次样条上采样 (Cubic Spline Upsampling)
% 对应文中公式 (4)
sig_upsampled = interp1(t_coarse, sig_discrete, t_fine, 'spline');

% 3. 计算互相关 (示意)
% 构造一个参考信号 (无时延)
ref_discrete = sig_fun(t_coarse, 0);
ref_upsampled = interp1(t_coarse, ref_discrete, t_fine, 'spline');

% 粗粒度互相关
[xc_coarse, lags_coarse] = xcorr(sig_discrete, ref_discrete);
lags_t_coarse = lags_coarse / fs_coarse;

% 高密度互相关 (对应文中公式 5)
[xc_fine, lags_fine] = xcorr(sig_upsampled, ref_upsampled);
lags_t_fine = lags_fine / fs_fine;

% 归一化幅度以便绘图
xc_coarse = xc_coarse / max(abs(xc_coarse));
xc_fine = xc_fine / max(abs(xc_fine));

%% === 2. JGR Format Plotting ===
figure('Units', 'inches', 'Position', [2, 2, 7, 6], 'Color', 'w');
t = tiledlayout(2, 1, 'TileSpacing', 'tight', 'Padding', 'compact');

% Global Style Settings
fontName = 'Arial';
fontSize = 11;
lineWidthMain = 1.5;
markerSize = 6;

% --- Panel (a): Signal Reconstruction ---
ax1 = nexttile; hold on; box on;
% 1. Plot True Analog Signal (Gray Dashed)
plot(t_fine*1e6, sig_true_continuous, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.2, ...
    'DisplayName', 'True Analog Waveform');
% 2. Plot Upsampled Reconstruction (Red Line)
plot(t_fine*1e6, sig_upsampled, 'r-', 'LineWidth', 1.2, ...
    'DisplayName', ['Spline Reconstruction (K=' num2str(K) ')']);
% 3. Plot Original Discrete Samples (Blue Stems/Dots)
% 使用 plot 画点，stem 有时在论文里显得太乱
plot(t_coarse*1e6, sig_discrete, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', markerSize, ...
    'DisplayName', 'Raw Discrete Samples');

% Zoom in on the peak area
xlim([0.35, 0.65]); 
ylabel('Amplitude (a.u.)', 'FontName', fontName, 'FontSize', fontSize);
set(gca, 'XTickLabel', []); % Hide X labels for top plot
grid on; set(gca, 'GridLineStyle', ':');

% Legend & Label
legend('Location', 'northeast', 'Box', 'off', 'FontName', fontName, 'FontSize', 9);
text(0.02, 0.92, '(a)', 'Units', 'normalized', 'FontName', fontName, 'FontSize', 12, 'FontWeight', 'bold');

% --- Panel (b): Cross-Correlation Improvement ---
ax2 = nexttile; hold on; box on;

% 1. Plot Fine Correlation (Red)
plot(lags_t_fine*1e6, xc_fine, 'r-', 'LineWidth', 1.5, ...
    'DisplayName', 'High-Density CCF');
% 2. Plot Coarse Correlation (Blue)
plot(lags_t_coarse*1e6, xc_coarse, 'bo--', 'MarkerFaceColor', 'b', 'MarkerSize', markerSize, ...
    'LineWidth', 1.0, 'DisplayName', 'Discrete CCF');

% Mark the True Peak vs Discrete Peak
[~, idx_max_fine] = max(xc_fine);
[~, idx_max_coarse] = max(xc_coarse);
plot(lags_t_fine(idx_max_fine)*1e6, xc_fine(idx_max_fine), 'k*', 'MarkerSize', 10, 'LineWidth', 1.5, ...
    'DisplayName', 'Global Max (Sub-sample)');

% Focus on the peak area
xlim([0.2, 0.5]); % Adjust based on the delay set above
ylim([0.8, 1.02]); % Focus on the peak top

xlabel('Time Lag (\mus)', 'FontName', fontName, 'FontSize', fontSize);
ylabel('Norm. Correlation', 'FontName', fontName, 'FontSize', fontSize);
grid on; set(gca, 'GridLineStyle', ':');

% Legend & Label
legend('Location', 'northeast', 'Box', 'off', 'FontName', fontName, 'FontSize', 9);
text(0.02, 0.92, '(b)', 'Units', 'normalized', 'FontName', fontName, 'FontSize', 12, 'FontWeight', 'bold');

% Global Font Adjustment
set(findall(gcf,'-property','FontSize'),'FontSize', fontSize);
set(findall(gcf,'-property','FontName'),'FontName', fontName);

fprintf('Figure 3 generated: Illustrating Spline Upsampling Principle.\n');