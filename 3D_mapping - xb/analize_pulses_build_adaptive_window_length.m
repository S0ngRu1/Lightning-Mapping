%% ========================================================================
%  脚本：分析脉冲密度并建立自适应窗口函数
% =========================================================================
% clear; clc; close all;
% 
% % --- 1. 加载和准备数据 ---
% % 加载您最长、最典型的滤波后信号
% % 您需要替换为您的真实数据读取方式
% signal_length = 2e8;
% r_loction = 3.65e8;
% waveform = filter_bp(read_signal('..\\20240822165932.6610CH1.dat', signal_length, r_loction), 30e6, 80e6, 5);
% noise_std = 3.545; % 假设这是已知的噪声水平
% threshold = noise_std * 3; % 假设用3倍标准差作为寻峰阈值
% 
% %% --- 2. 步骤一：执行数据普查 (滑动窗口扫描) ---
% fprintf('开始执行脉冲密度普查...\n');
% scout_window_len = 20000; % 侦察窗口大小
% scout_step = 5000;       % 侦察窗口的滑动步长
% min_peak_dist_in_scan = 200; % 侦察时使用的大致脉冲最小间距
% 
% % 计算滑动窗口的数量
% num_windows = floor((length(waveform) - scout_window_len) / scout_step) + 1;
% 
% % 初始化结果容器
% window_centers = zeros(num_windows, 1);
% pulse_counts = zeros(num_windows, 1);
% 
% h = waitbar(0, '正在进行密度扫描...');
% for i = 1:num_windows
%     waitbar(i/num_windows, h);
%     
%     win_start = (i-1)*scout_step + 1;
%     win_end = win_start + scout_window_len - 1;
%     
%     segment = waveform(win_start:win_end);
%     
%     [~, locs] = findpeaks(abs(segment), 'MinPeakHeight', threshold, 'MinPeakDistance', min_peak_dist_in_scan);
%     
%     window_centers(i) = win_start + scout_window_len/2;
%     pulse_counts(i) = numel(locs);
% end
% close(h);
% fprintf('密度普查完成。\n');
% 
% %% --- 3. 步骤二：数据可视化分析 ---
% 
% % 图1：密度随时间的变化图
% figure('Name', '脉冲密度分析');
% subplot(2,1,1);
% plot(window_centers+r_loction, pulse_counts, '.-');
% title('信号脉冲密度随时间变化图');
% xlabel('采样点位置');
% ylabel(sprintf('在 %d 点窗口内的脉冲数', scout_window_len));
% grid on;
% 
% % 图2：密度分布直方图 (最关键的图)
% subplot(2,1,2);
% histogram(pulse_counts, 'BinMethod', 'integers');
% title('脉冲密度分布直方图');
% xlabel(sprintf('在 %d 点窗口内的脉冲数', scout_window_len));
% ylabel('频次');
% grid on;
% set(gca, 'YScale', 'log'); % 使用对数Y轴可能看得更清楚
% 
% 
% 
% 


% ------------------------------------------------------------------------

%% ========================================================================
%  脚本：使用 find_pulses_advanced 进行高精度脉冲密度普查
% =========================================================================
clear; clc; close all;

% --- 1. 加载和准备数据 ---
signal_length = 2e8;
r_loction = 3.65e8;
sampling_rate_hz = 200e6; % 必须定义采样率
fprintf('正在读取和滤波信号 (长度: %d)...\n', signal_length);
waveform = filter_bp(read_signal('..\\20240822165932.6610CH1.dat', signal_length, r_loction), 30e6, 80e6, 5);
noise_std = 3.545; % 已知的噪声水平
detection_threshold_factor = 3; % 检测因子

%% --- 2. 步骤一：执行高精度密度普查 (滑动窗口扫描) ---
fprintf('开始执行高精度脉冲密度普查 (这将非常耗时)...\n');
scout_window_len = 20000; % 侦察窗口大小
scout_step = 5000;       % 侦察窗口的滑动步长

% 计算滑动窗口的数量
num_windows = floor((length(waveform) - scout_window_len) / scout_step) + 1;

% 初始化结果容器
window_centers = zeros(num_windows, 1);
pulse_counts = zeros(num_windows, 1);

h = waitbar(0, '正在进行高精度密度扫描 (可能需要数小时)...');
for i = 1:num_windows
    waitbar(i/num_windows, h, sprintf('处理窗口 %d/%d...', i, num_windows));
    
    win_start = (i-1)*scout_step + 1;
    win_end = win_start + scout_window_len - 1;
    segment = waveform(win_start:win_end);
    pulse_catalog_segment = find_pulses_advanced(segment, noise_std, sampling_rate_hz, detection_threshold_factor);
    pulse_counts(i) = numel(pulse_catalog_segment);
    window_centers(i) = win_start + scout_window_len/2;
end
close(h);
fprintf('高精度密度普查完成。\n');

%% --- 3. 步骤二：数据可视化分析 ---
% 图1：密度随时间的变化图
figure('Name', '高精度脉冲密度分析');
subplot(2,1,1);
plot(window_centers + r_loction, pulse_counts, '.-');
title('高精度信号脉冲密度随时间变化图');
xlabel('采样点绝对位置');
ylabel(sprintf('在 %d 点窗口内的精确脉冲数', scout_window_len));
grid on;
% 图2：密度分布直方图
subplot(2,1,2);
histogram(pulse_counts, 'BinMethod', 'integers');
title('高精度脉冲密度分布直方图');
xlabel(sprintf('在 %d 点窗口内的精确脉冲数', scout_window_len));
ylabel('频次');
grid on;
set(gca, 'YScale', 'log'); % 使用对数Y轴可能看得更清楚


fprintf('开始执行全自动区域分析与四档建模...\n');

%% --- 2. 步骤一 & 二：周期性扫描，自动识别最稀疏和最密集区域 ---

analysis_chunk_size = 2e7; % 按2千万个采样点为一块进行分析
num_chunks = floor(length(waveform) / analysis_chunk_size);

% 初始化用于存储每个区块分析结果的变量
chunk_avg_counts = zeros(num_chunks, 1);
chunk_start_locs = zeros(num_chunks, 1);

fprintf('正在以 %d 点为单位，分块扫描信号...\n', analysis_chunk_size);
for i = 1:num_chunks
    chunk_start = (i-1) * analysis_chunk_size + 1;
    chunk_end = i * analysis_chunk_size;
    
    % 找到落在当前区块内的普查窗口的索引
    indices_in_chunk = find(window_centers >= chunk_start & window_centers < chunk_end);
    
    if ~isempty(indices_in_chunk)
        % 计算该区块的平均脉冲数
        chunk_avg_counts(i) = mean(pulse_counts(indices_in_chunk));
    else
        chunk_avg_counts(i) = 0; % 如果区块内没有普查点，则密度为0
    end
    chunk_start_locs(i) = chunk_start + r_loction;
end

% 找到平均脉冲数最大（最密集）的区块
[max_avg_count, max_idx] = max(chunk_avg_counts);
densest_zone_loc = chunk_start_locs(max_idx);

% 找到平均脉冲数最小（最稀疏）的区块，要排除全为0的区块
non_zero_chunks_indices = find(chunk_avg_counts > 0);
[min_avg_count, min_idx_relative] = min(chunk_avg_counts(non_zero_chunks_indices));
min_idx_absolute = non_zero_chunks_indices(min_idx_relative);
sparsest_zone_loc = chunk_start_locs(min_idx_absolute);

fprintf('\n--- 自动区域分析结果 ---\n');
fprintf('最密集区域位于 %.2e 附近，平均脉冲数为: %.2f\n', densest_zone_loc, max_avg_count);
fprintf('最稀疏区域位于 %.2e 附近，平均脉冲数为: %.2f\n', sparsest_zone_loc, min_avg_count);


%% --- 3. 步骤三：基于分析结果，自动划分四档阈值 ---

% 我们需要确定三个阈值 T1, T2, T3 来划分四个区间

T1 = round(min_avg_count );
T3 = round(max_avg_count );

% 确保 T1 < T3
if T1 >= T3
    fprintf('警告: 密集区与稀疏区密度差异过小，无法清晰划分四档，将采用三档划分。\n');
    T2 = T1; % 让T1和T2相等，实际只有三档
else
    % T2可以取T1和T3的对数中心点，这在处理跨度大的数据时更合理
    T2 = round(10^( (log10(T1) + log10(T3)) / 2 ));
end

fprintf('\n--- 自动分档阈值确定 ---\n');
fprintf('级别1 (非常稀疏, 窗口 4096): 脉冲数 <= %d\n', T1);
fprintf('级别2 (较为稀疏, 窗口 2048): 脉冲数 > %d 且 <= %d\n', T1, T2);
fprintf('级别3 (较为密集, 窗口 1024): 脉冲数 > %d 且 <= %d\n', T2, T3);
fprintf('级别4 (非常密集, 窗口 512): 脉冲数 > %d\n\n', T3);

