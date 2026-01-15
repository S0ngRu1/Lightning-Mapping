%% ================== 0. 初始化与参数设置 ==================
clear; clc; close all;

% --- 核心参数设置 ---
r_location_yld = 3.708e8;       % Start location (Sample index)
signal_length = 1e6;            % Read length (Sample points, 5e5 = 2.5ms)
sampling_interval_ns = 5;       % Sampling interval 5ns
fs = 1e9 / sampling_interval_ns;% Sampling rate 200MHz


% --- 筛选阈值参数 (来自您的要求) ---
Error_Threshold = 1;            % 误差阈值
Rcorr_Default = 0.3;            % 默认相关系数阈值

% --- 文件路径 (请根据实际情况确认) ---
file_ch1 = '..\20240822165932.6610CH1.dat';
file_ch4 = '..\20240822165932.6610CH4.dat';
file_result = 'results\result_yld_3.65e8_5.6e8_window_ADAPTIVE_1e4_factor3_with_error.txt';

% --- 时间轴构建 ---
time_conversion_factor = sampling_interval_ns * 1e-3; % 转换为微秒
x_indices = 0 : signal_length - 1;
time_us = x_indices * time_conversion_factor;

%% ================== 1. 信号读取与处理 ==================
fprintf('正在读取和处理波形信号...\n');

% 1. 读取原始数据
ch1_raw = read_signal(file_ch1, signal_length, r_location_yld);
ch4_raw = read_signal(file_ch4, signal_length, r_location_yld);

% 2. 通道1 (VHF) 处理: 带通滤波 30-80MHz
ch1_proc = filter_bp(ch1_raw, 30e6, 80e6, fs, 5);

% 3. 通道4 (快电场) 处理: 基线去除 + 小波去噪
% 3. 通道4 (快电场) 处理: 组合降噪
% 3.1 基线去除 (移动中值)
baseline = movmedian(ch4_raw, 1024);
E_temp = double(ch4_raw - baseline);

% 3.2 第一步：低通滤波 (去除高频底噪) -> 最有效的一步
% 截止频率建议：5MHz ~ 10MHz。如果脉冲变宽太严重，请调高此值。
f_lp = 8e6; % 8 MHz
[b_lp, a_lp] = butter(4, f_lp / (fs/2), 'low');
E_lp = filtfilt(b_lp, a_lp, E_temp);

% 3.3 第二步：Savitzky-Golay 平滑 (去除残留毛刺，保留峰值)
% 窗口长度 21 比较适中
E_fast_denoised = sgolayfilt(E_lp, 3, 21);

% (可选) 如果幅值太小看不清，可以乘以一个增益系数
% E_fast_denoised = E_fast_denoised * 2;

%% ================== 2. 辐射源数据加载与严格筛选 ==================
fprintf('正在读取和筛选定位数据...\n');
if ~isfile(file_result), error('文件不存在: %s', file_result); end
opts = detectImportOptions(file_result); opts.VariableNamingRule = 'preserve';
T = readtable(file_result, opts);

% --- 2.1 几何与物理范围筛选 ---
% 仅保留当前绘图时间窗口内的数据
filter_loc_min = r_location_yld;
filter_loc_max = r_location_yld + signal_length;

valid_geo = T.Start_loc >= filter_loc_min & ...
            T.Start_loc <= filter_loc_max & ...
            T.Elevation < 18 & ...
            T.Elevation > 10 & ...
            abs(T.t123) < 1;

% --- 2.2 坏点区域剔除 (根据您的代码逻辑) ---
bad_region = false(height(T), 1);
if ismember('Azimuth', T.Properties.VariableNames)
    % 注意：您提供的条件 (El>0 & El<0) 实际上不会剔除任何点，这里保留代码逻辑
    bad_region = (T.Azimuth > 160) & (T.Azimuth < 250) & (T.Elevation > 0) & (T.Elevation < 0);
end

% --- 2.3 自适应相关系数筛选 (核心逻辑修正) ---
mask_rcorr = false(height(T), 1);
if ismember('Win_len', T.Properties.VariableNames)
    idx_512 = (T.Win_len == 512);   mask_rcorr(idx_512) = T.Rcorr(idx_512) > 0.5;
    idx_1024 = (T.Win_len == 1024); mask_rcorr(idx_1024) = T.Rcorr(idx_1024) > 0.3;
    idx_2048 = (T.Win_len == 2048); mask_rcorr(idx_2048) = T.Rcorr(idx_2048) > 0.2;
    idx_4096 = (T.Win_len == 4096); mask_rcorr(idx_4096) = T.Rcorr(idx_4096) > 0.1;
    other = ~ismember(T.Win_len, [512, 1024, 2048, 4096]);
    mask_rcorr(other) = T.Rcorr(other) > Rcorr_Default;
else
    mask_rcorr = T.Rcorr > Rcorr_Default;
end

% --- 2.4 误差阈值筛选 ---
if ismember('Err_Az', T.Properties.VariableNames)
    mask_error = (T.Err_Az < Error_Threshold) & (T.Err_El < Error_Threshold);
else
    mask_error = true(height(T), 1);
end

% --- 2.5 应用综合筛选 ---
final_idx = valid_geo & (~bad_region) & mask_rcorr & mask_error;
filteredTable = T(final_idx, :);

% 计算相对时间 (us) 和 仰角
vhf_time = (filteredTable.Start_loc - r_location_yld) * time_conversion_factor;
vhf_el = filteredTable.Elevation;
vhf_color = vhf_time; % 颜色随时间变化

fprintf('筛选完成，保留点数: %d\n', height(filteredTable));

%% ================== 3. JGR 风格绘图 ==================
% 设置画布 (16cm x 14cm)
figure('Units', 'centimeters', 'Position', [5, 5, 16, 14], 'Color', 'w');
t = tiledlayout(2, 1, 'TileSpacing', 'tight', 'Padding', 'compact');

% 字体设置
font_name = 'Arial';
font_size = 9;

% --- 子图 1: 信号波形对比 (双Y轴) ---
ax1 = nexttile;
% 左轴: VHF (CH1)
yyaxis(ax1, 'left');
plot(ax1, time_us, ch1_proc, 'Color', [0, 0.447, 0.741], 'LineWidth', 0.6); % 蓝色
ylabel('VHF Signal (ADU)', 'FontName', font_name, 'FontSize', font_size);
set(ax1, 'YColor', [0, 0.447, 0.741]); 

% 右轴: 快电场 (CH4)
yyaxis(ax1, 'right');
plot(ax1, time_us, E_fast_denoised, 'Color', [0.85, 0.325, 0.098], 'LineWidth', 0.6); % 红色
ylabel('Fast E-field (ADU)', 'FontName', font_name, 'FontSize', font_size);
set(ax1, 'YColor', [0.85, 0.325, 0.098]);

title('(a) VHF and Fast Electric Field Waveforms', 'FontName', font_name, 'FontSize', font_size+1, 'FontWeight', 'bold');
grid(ax1, 'on');
apply_jgr_style(ax1, font_name, font_size);
set(ax1, 'XTickLabel', []); % 隐藏X轴标签以紧凑布局
xlim(ax1, [0, max(time_us)]);

% --- 子图 2: 仰角随时间变化 ---
ax2 = nexttile;
point_size = 15;
% 绘制散点图
scatter(ax2, vhf_time, vhf_el, point_size, vhf_color, 'filled', 'MarkerFaceAlpha', 0.7);

ylabel('Elevation (°)', 'FontName', font_name, 'FontSize', font_size);
xlabel('Time (\mus)', 'FontName', font_name, 'FontSize', font_size);
title('(b) Evolution of VHF Source Elevation', 'FontName', font_name, 'FontSize', font_size+1, 'FontWeight', 'bold');

% Colorbar 设置
cb = colorbar;
cb.Label.String = 'Time (\mus)';
cb.Label.FontSize = font_size;
cb.Label.FontName = font_name;
colormap(ax2, 'jet');

grid(ax2, 'on');
apply_jgr_style(ax2, font_name, font_size);
xlim(ax2, [0, max(time_us)]);

% 链接X轴 (缩放同步)
linkaxes([ax1, ax2], 'x');

fprintf('绘图完成。\n');

%% ================== 辅助函数定义 ==================

function apply_jgr_style(ax, fname, fsize)
    % 统一应用 JGR 期刊风格
    set(ax, 'FontName', fname, 'FontSize', fsize);
    set(ax, 'LineWidth', 1.0);      % 坐标轴线宽
    set(ax, 'TickDir', 'in');       % 刻度朝内
    set(ax, 'Box', 'on');           % 开启边框
    set(ax, 'GridAlpha', 0.15);     % 网格透明度
    set(ax, 'GridLineStyle', '--'); % 网格样式
    set(ax, 'Layer', 'top');        % 刻度浮于数据之上
end

function signal = read_signal(signal_path, r_length, r_location)
    % 读取二进制信号文件
    fid = fopen(signal_path, 'r');
    if fid == -1, error('无法打开文件: %s', signal_path); end
    % r_location 为采样点位置，int16占2字节
    status = fseek(fid, r_location * 2, 'bof'); 
    if status == -1, fclose(fid); error('fseek 失败'); end
    signal = fread(fid, r_length, 'int16');
    fclose(fid);
end

function filtered_signal = filter_bp(signal, f1, f2, fs, order)
    % 巴特沃斯带通滤波
    fn = fs/2;
    Wn = [f1 f2]/fn;
    [b,a] = butter(order, Wn);
    filtered_signal = filtfilt(b, a, double(signal));
end