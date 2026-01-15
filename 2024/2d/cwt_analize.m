%% === 参数设置 ===
fs = 200e6;
signal_length = 100;
%正先导信号
% r_loction_yld = 383000000 + 1400000 - 30 + 150 + 4848;

%负先导信号
r_loction_yld = 367000000+2800000+760;

%% === 1. 数据读取与处理 (仅通道1) ===
% 指定只读取 CH1
target_channel = 'CH1';
filename = sprintf('..\\20240822165932.6610%s.dat', target_channel);

% 读取信号 (确保 read_signal 函数在路径中)
signal = read_signal(filename, signal_length, r_loction_yld);

% 滤波处理 (确保 filter_bp 函数在路径中)
signal = filter_bp(signal, 30e6, 80e6, 5);

% 构建滤波器组
fb = cwtfilterbank( ...
    'SignalLength', length(signal), ...
    'SamplingFrequency', fs, ...
    'Wavelet', 'morse', ...
    'FrequencyLimits', [25e6 85e6], ...   
    'VoicesPerOctave', 48);              

% 进行小波变换
[cfs, f] = wt(fb, signal);

%% === 2. 绘图展示 ===
figure('Name', 'CH1 信号分析', 'Color', 'w', 'Position', [500, 200, 600, 500]);

% 准备时间轴数据 (us)
t_us = (1:length(signal)) / fs * 1e6;

% --- 子图 1: 信号时域波形图 (原顺序互换) ---
subplot(2, 1, 1);
plot(t_us, signal, 'b', 'LineWidth', 1);  
xlabel('时间 (\mus)');
ylabel('信号幅度');
title('信号时域波形图'); % 按要求修改标题
grid on;
xlim([min(t_us), max(t_us)]); % 保持X轴紧凑

% --- 子图 2: 信号时频图(CWT) ---
subplot(2, 1, 2);
scalogram_data = abs(cfs);
imagesc(t_us, f/1e6, scalogram_data);  
axis xy; % 频率轴方向修正：低频在下，高频在上
colormap('jet'); % 使用 jet 色图使能量分布更清晰
xlabel('时间 (\mus)');
ylabel('频率 (MHz)');
title('信号时频图(CWT)'); % 按要求修改标题
colorbar;

%% === 辅助函数占位符 (请确保你的脚本文件夹里包含这两个函数) ===
% 如果你的脚本是单个文件运行，请保留之前的 read_signal 和 filter_bp 函数定义
% function signal = read_signal(...)
% function filtered_signal = filter_bp(...)