clc;
clear;
close all;
% --- 数据文件路径 ---
signal_file = '..\20240822165932.6610CH4.dat';
vhf_file = 'results\20240822165932_result_yld_3.6e8_5.6e8_window_1024_256_阈值4倍标准差_去零飘_30_80_hann.txt';
% --- 加载位置参数 ---
r_loction_yld = 3.83e8;    % 加载起始点
signal_length = 0.4e7;     % 加载长度 

% --- 梯级检测参数  ---
MIN_PROMINENCE = 0.1;   % 最小突起度 (0-1之间)，值越大检测越严，抗噪越好
MIN_DISTANCE_MS = 0.005; % 最小梯级间距(ms)，5us防止检测到同一次脉冲的振铃
BASELINE_WINDOW = 1024;   % 去趋势的滑动窗口大小

% --- 绘图参数 ---
POINT_SIZE = 15;         % VHF散点大小


%% ================== 1. 数据加载与预处理 ==================
disp('1. 正在加载和预处理数据...');

% --- 1a. 加载快电场信号 ---
if exist('read_signal', 'file')
    ch1_yld = read_signal(signal_file, signal_length, r_loction_yld);
else
    error('未找到 read_signal 函数，请确保该函数在当前工作目录或路径中。');
end

% 信号归一化
y_abs_max = max(abs(ch1_yld));
if y_abs_max ~= 0
    ch1_normalized = ch1_yld / y_abs_max;
else
    ch1_normalized = ch1_yld;
end

% 创建时间轴 (ms)
sampling_interval_ns = 5;
time_conversion_factor = sampling_interval_ns * 1e-6; % ns 转 ms
time_ms = (0 : signal_length - 1) * time_conversion_factor;

% --- 1b. 加载 VHF 辐射源数据 ---
result1 = readtable(vhf_file);
% 过滤数据
filter_loc_min = r_loction_yld;
filter_loc_max = r_loction_yld + signal_length;
logicalIndex = abs(result1.t123) < 1 & abs(result1.Rcorr) > 0.6 & ...
    result1.Start_loc < filter_loc_max & result1.Start_loc > filter_loc_min;
filteredTable1 = result1(logicalIndex, :);

% 转换 VHF 时间
vhf_time_ms = (filteredTable1.Start_loc - r_loction_yld) * time_conversion_factor;
vhf_elevation = filteredTable1.Elevation;


%% ================== 2. 闪电先导梯级参数求解 ==================
disp('2. 正在执行梯级脉冲检测...');

% 2a. 信号去趋势 
baseline = movmedian(ch1_normalized, BASELINE_WINDOW);
E_fast = ch1_normalized - baseline;
% 2b. 脉冲检测
E_to_detect = abs(E_fast);

[pks, locs_time, widths_ms, proms] = findpeaks(E_to_detect, time_ms, ...
    'MinPeakProminence', MIN_PROMINENCE, ...
    'MinPeakDistance', MIN_DISTANCE_MS, ...
    'WidthReference', 'halfheight');

% 2c. 计算物理参数
leader_step_times = locs_time;                % 梯级发生时刻 (ms)
step_durations_us = widths_ms * 1000;         % 持续时间 (us)
step_intervals_us = diff(locs_time) * 1000;   % 间隔时间 (us)

% 2d. 打印统计结果
% 2d. 打印统计结果
% 2d. 打印统计结果
fprintf('--------------------------------------------------------\n');
fprintf('检测参数         : Prominence=%.2f, MinDist=%.3f ms\n', MIN_PROMINENCE, MIN_DISTANCE_MS);
fprintf('梯级检测统计结果\n');
% 打印信号段范围 (使用 %.0f 以显示完整的整数索引，避免科学计数法丢失精度)
fprintf('信号段范围: %.4e - %.4e (信号长度: %.1f ms)\n', r_loction_yld, r_loction_yld + signal_length, signal_length*5/1e6);
fprintf('--------------------------------------------------\n');
fprintf('检测到梯级总数   : %d 个\n', length(locs_time));
if ~isempty(step_durations_us)
    fprintf('梯级持续时间 (FWHM): 平均=%.2f us, 中位数=%.2f us, 标准差=%.2f us, [最小=%.2f, 最大=%.2f] us\n', ...
        mean(step_durations_us), median(step_durations_us), std(step_durations_us), min(step_durations_us), max(step_durations_us));
    fprintf('梯级间隔时间       : 平均=%.2f us, 中位数=%.2f us, 标准差=%.2f us, [最小=%.2f, 最大=%.2f] us\n', ...
        mean(step_intervals_us), median(step_intervals_us), std(step_intervals_us), min(step_intervals_us), max(step_intervals_us));
else
    fprintf('未检测到足够的梯级，请尝试减小 MIN_PROMINENCE 阈值。\n');
end
fprintf('--------------------------------------------------------\n');

%% ================== 3. 绘图 A: 主图 (融合图 + 叠加检测结果) ==================
fig1 = figure('Color', [1 1 1], 'Position', [50, 100, 1000, 600], 'Name', 'VHF与电场综合分析');
ax = gca;

% --- 左轴: 快电场 ---
yyaxis(ax, 'left');
% 绘制原始归一化电场 (黑色细线)
hE = plot(ax, time_ms, ch1_normalized, 'Color', 'k', 'LineWidth', 1.2); hold on;
% [新增] 叠加红点，显示算法检测到的梯级位置 (用于人工核对)
% 为了让红点显示在原始波形上，我们需要找到 locs_time 对应的原始信号幅值
% 使用 interp1 插值找到对应时间点的原始电场值
pks_on_raw_signal = interp1(time_ms, ch1_normalized, locs_time, 'nearest');
hDet = plot(ax, locs_time, pks_on_raw_signal, 'ro', 'MarkerSize', 6, 'LineWidth', 1.5, 'DisplayName', '自动检测的梯级');

ylabel(ax, '归一化电场强度', 'FontSize', 12, 'FontWeight', 'bold');
set(ax, 'YColor', 'k'); % 左轴黑色
ylim(ax, [-1.1, 1.1]);

% --- 右轴: VHF仰角 ---
yyaxis(ax, 'right');
scatter(ax, vhf_time_ms, vhf_elevation, ...
    POINT_SIZE, vhf_time_ms, 'filled', 'MarkerFaceAlpha', 0.7);
ylabel(ax, '仰角 (Elevation / °)', 'FontSize', 12, 'FontWeight', 'bold');
set(ax, 'YColor', [0.2 0.2 0.2]); % 右轴深灰色
ylim(ax, [0, 85]); yticks(ax, 0:10:85);

% --- 公共轴与修饰 ---
xlabel(ax, '时间 (ms)', 'FontSize', 12, 'FontWeight', 'bold');
time_min_ms = 0;
time_max_ms = (signal_length - 1) * time_conversion_factor;
xlim(ax, [time_min_ms, time_max_ms]);

colormap(ax, 'parula');
hcb = colorbar(ax);
ylabel(hcb, 'VHF辐射时间 (ms)', 'FontSize', 11);
clim(ax, [time_min_ms, time_max_ms]);

title('闪电快电场梯级检测与VHF辐射源时序图', 'FontSize', 16, 'FontWeight', 'bold');
grid(ax, 'on');
set(ax, 'FontSize', 11, 'LineWidth', 1.2, 'Box', 'on', 'GridAlpha', 0.3, 'GridLineStyle', '--');

% 添加图例 (仅显示左轴的电场和检测点)
legend([hE, hDet], {'快电场波形', '检测到的梯级(Steps)'}, 'Location', 'northwest', 'FontSize', 10);


%% ================== 4. 绘图 B: 梯级参数统计图 ==================
if ~isempty(step_durations_us)
    fig2 = figure('Color', [1 1 1], 'Position', [100, 100, 800, 400], 'Name', '梯级参数统计');

    % 间隔时间直方图
    subplot(1, 2, 1);
    histogram(step_intervals_us, 30, 'FaceColor', [0 0.4470 0.7410]);
    xlabel('梯级间隔时间 (\mus)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('频次 (Count)', 'FontSize', 12, 'FontWeight', 'bold');
    title(['平均间隔: ' sprintf('%.2f', mean(step_intervals_us)) ' \mus'], 'FontSize', 14);
    grid on;

    % 持续时间直方图
    subplot(1, 2, 2);
    histogram(step_durations_us, 30, 'FaceColor', [0.8500 0.3250 0.0980]);
    xlabel('梯级持续时间 (FWHM, \mus)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('频次 (Count)', 'FontSize', 12, 'FontWeight', 'bold');
    title(['平均持续: ' sprintf('%.2f', mean(step_durations_us)) ' \mus'], 'FontSize', 14);
    grid on;
end