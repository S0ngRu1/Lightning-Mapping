%% ========================================================================
%  脚本：对真实闪电信号进行小波变换动态频谱分析
% =========================================================================
clear; clc; close all;

signal_length = 8e6;        
r_loction_yld = 4.69e8;     % 起始位置
fs = 200e6;                 % 采样率 200 MHz
ts = 1/fs;                  % 采样周期

ch1_yld = read_signal('..\\20240822165932.6610CH1.dat', signal_length, r_loction_yld);
filtered_signal = filter_bp(ch1_yld, 30e6, 80e6, 5);


%% --- 3. 使用 cwt 进行动态频谱分析 ---

% 创建时间轴，单位：微秒 (μs)
time_us = (0:signal_length-1) * ts * 1e6; 

% 手动计算CWT
[cfs, f] = cwt(filtered_signal, fs, ...
                'Wavelet', 'morse', ...                 % 使用灵活的 Morse 小波
                'FrequencyLimits', [30e6, 80e6], ...      
                'VoicesPerOctave', 16);                % 提高一点频率分辨率


figure('Name', '闪电信号动态频谱 - 微秒时间轴', 'Position', [150, 150, 1200, 600]);

% 使用 pcolor 手动绘图
pcolor(time_us, f/1e6, abs(cfs)); % 将频率转换为 MHz
shading interp;

% 美化图形
title(sprintf('闪电信号动态频谱图 (Scalogram) @ %.2e', r_loction_yld), 'FontSize', 14);
xlabel('时间 (微秒, μs)', 'FontSize', 12);
ylabel('频率 (MHz)', 'FontSize', 12);
set(gca, 'YScale', 'log'); % 对数频率轴可以更好地展示细节
colormap('jet');
h_bar = colorbar;
ylabel(h_bar, '能量幅度');