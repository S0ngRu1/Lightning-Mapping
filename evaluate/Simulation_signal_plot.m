%% 仿真信号生成流程可视化脚本
clear; clc; close all;

%% === 1. 环境与参数初始化 (为您的代码片段提供必要的上下文) ===
% 为了让您的代码运行，我们需要定义以下变量：
% N, sig_len, snr_val, b, a, win_len, phase_shift

fs = 200e6;                     % 采样率 200 MHz (用于时间轴显示)
N = 4096;                       % 总仿真点数 (背景噪声长度)
sig_len = 512;                  % 有效信号长度 (Burst长度)
win_len = 1024;                 % 最终截取的分析窗口长度
snr_val = 15;                   % 信噪比 (dB)，设置高一点以便观察
true_delay_samples = 20.5;      % 设置一个明显的亚采样延迟 (20.5个点)

% 设计一个带通滤波器 (例如 30MHz - 80MHz)
f_low = 30e6; f_high = 80e6;
[b, a] = butter(4, [f_low, f_high]/(fs/2));

% 预计算频域相移向量 (用于实现亚采样延迟)
% 使用双边频率轴以确保 IFFT 的结果是实数
f_axis = [0:ceil(N/2)-1, -floor(N/2):-1] * (fs/N);
phase_shift = exp(-1j * 2 * pi * f_axis * (true_delay_samples/fs));

fprintf('参数准备完毕，开始执行您的仿真代码片段...\n');

%% === 2. 核心代码片段 (直接嵌入您提供的内容) ===

% 1. 信号生成
raw_sig = randn(1, N);
mask = zeros(1, N);
center = floor(N/2);
start_p = max(1, center - floor(sig_len/2));
end_p = min(N, start_p + sig_len - 1);
mask(start_p:end_p) = 1;
sig = raw_sig .* mask;

% 2. 延迟与噪声
SIG = fft(sig);
% 注意：为了让 'symmetric' 正常工作，输入需要是共轭对称的，这由上面的双边 f_axis 保证
sig_delayed = ifft(SIG .* phase_shift, 'symmetric');

% 计算缩放因子 (基于单位噪声方差和目标SNR)
scale = sqrt(1 * 10^(snr_val/10));

r1 = sig * scale + randn(1, N);
r2 = sig_delayed * scale + randn(1, N);

% 3. 滤波与截窗
% 为了可视化滤波效果，我们先把滤波后的全长信号保存下来
r1_filt_full = filter(b, a, r1);
r2_filt_full = filter(b, a, r2);

w_start = floor(N/2) - floor(win_len/2);
idx = w_start : (w_start + win_len - 1);
% 确保索引不越界
idx = idx(idx >= 1 & idx <= N);

% 最终截取的窗口信号 (您的代码中的 w1, w2)
w1 = r1_filt_full(idx);
w2 = r2_filt_full(idx);


%% === 3. 信号可视化绘图 ===
fprintf('仿真完成，正在绘制结果...\n');

% 创建时间轴用于显示 (微秒)
t_full = (0:N-1) / fs * 1e6;
t_win = (idx-1) / fs * 1e6;

figure('Units', 'pixels', 'Position', [100, 100, 1200, 900], 'Color', 'w');
tiledlayout(3, 1, 'TileSpacing', 'tight', 'Padding', 'compact');

% --- 子图1：原始源信号生成 (Step 1) ---
nexttile; hold on; grid on;
plot(t_full, raw_sig, 'Color', [0.8 0.8 0.8], 'DisplayName', '背景噪声 (raw\_sig)');
plot(t_full, sig, 'b', 'LineWidth', 1.5, 'DisplayName', '生成的源信号 Burst (sig)');
plot(t_full, mask * max(abs(sig))*1.2, 'r--', 'LineWidth', 1, 'DisplayName', '掩码 (mask)');
ylabel('幅值'); title('源信号生成');
legend; xlim([min(t_full), max(t_full)]);
ylim([-4, 5]); % 限制一下Y轴看清中间

% --- 子图2：接收到的含噪信号 (Step 2) ---
nexttile; hold on; grid on;
% 为了看清细节，只放大显示中间一部分区域
zoom_range = center-win_len : center+win_len; 
zoom_range = zoom_range(zoom_range > 0 & zoom_range <= N);

plot(t_full(zoom_range), r1(zoom_range), 'Color', [0.6 0.8 1], 'LineWidth', 0.5, 'DisplayName', '接收信号 r1 (含噪)');
% 为了凸显延迟，把 r2 画在下面一点
plot(t_full(zoom_range), r2(zoom_range) - 10, 'Color', [1 0.6 0.6], 'LineWidth', 0.5, 'DisplayName', '接收信号 r2 (延迟+含噪，偏移显示)');

% 叠加缩放后的纯净信号作为参考
plot(t_full(zoom_range), sig(zoom_range)*scale, 'b--', 'LineWidth', 1.5, 'DisplayName', '纯净参考 r1');
plot(t_full(zoom_range), sig_delayed(zoom_range)*scale - 10, 'r--', 'LineWidth', 1.5, 'DisplayName', '纯净参考 r2 (延迟)');

ylabel('幅值 (偏移显示)'); 
title(['含噪信号 (SNR=' num2str(snr_val) 'dB, 延迟=' num2str(true_delay_samples) '样点)']);
legend; xlim([t_full(zoom_range(1)), t_full(zoom_range(end))]);

% --- 子图3：最终滤波和截窗后的信号 (Step 3) ---
nexttile; hold on; grid on;
plot(t_win, w1, 'b', 'LineWidth', 1.5, 'DisplayName', '最终窗口信号 w1 (参考)');
plot(t_win, w2, 'r', 'LineWidth', 1.0, 'DisplayName', '最终窗口信号 w2 (延迟)');

xlabel('时间 (\mus)'); ylabel('幅值'); 
title(['滤波后的信号 (窗口长度=' num2str(win_len) ')']);
legend; xlim([min(t_win), max(t_win)]);

sgtitle('仿真信号生成流程可视化');