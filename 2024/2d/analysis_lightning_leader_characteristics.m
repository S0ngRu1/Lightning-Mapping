clear; 
close all; 
clc;

%% ==================== 1. 参数设置与数据加载 ====================
% --- 用户需根据实际情况修改的参数 ---
DATA_FILE = '20230718175104_result_yld_3e8_6e8_window_512_128_阈值4倍标准差_去零飘1_30_80_hann.txt'; % 您的数据文件名
SAMPLING_RATE = 200e6;            % 数据采集卡采样率 (Hz), 200 MS/s
ASSUMED_HEIGHT = 10;            % 假设的放电平均高度 (米), 用于将角度转换成距离

% --- 速度计算参数 ---
TIME_WINDOW_VEL = 0.00010;        % 速度计算的时间窗口 (秒)

% --- 列定义  ---
COL_START_LOC = 1;
COL_AZIMUTH   = 8;
COL_ELEVATION = 9;
COL_RCORR     = 10;
COL_T123      = 11;

% --- 定义两个事件的时间窗口 (单位: 采样点) ---
event1_range = [5.2e8, 5.3e8]; % 用于左子图
event2_range = [4.0e8, 4.1e8]; % 用于右子图

% --- 加载原始数据 ---
fprintf('正在加载原始数据: %s\n', DATA_FILE);
try
    data_raw = readmatrix(DATA_FILE);
catch
    error('数据文件加载失败，请检查文件名或文件格式是否正确。');
end
fprintf('原始数据点数: %d\n', size(data_raw, 1));


%% ==================== 2. 分别处理两个事件 ====================

% --- 处理事件1 (右子图数据) ---
fprintf('\n--- 正在处理事件1 (%.2e to %.2e) ---\n', event1_range(1), event1_range(2));
velocities1 = []; time_midpoints1 = []; avg_vel1 = NaN; % 初始化结果
% 筛选事件1数据
logicalIndex1 = ...
    abs(data_raw(:, COL_RCORR)) > 0.6 & ...
    data_raw(:, COL_START_LOC) > event1_range(1) & ...
    data_raw(:, COL_START_LOC) < event1_range(2) & ...
    data_raw(:, COL_ELEVATION) < 80 & ...
    abs(data_raw(:, COL_T123)) < 1;
data1 = data_raw(logicalIndex1, :);
fprintf('事件1筛选后有效数据点数: %d\n', size(data1, 1));

if ~isempty(data1)
    % 预处理事件1数据
    time_samples1 = data1(:, COL_START_LOC);
    azimuth_deg1 = data1(:, COL_AZIMUTH);   % 【修正】使用常量
    elevation_deg1 = data1(:, COL_ELEVATION); % 【修正】使用常量
    time_sec1 = (time_samples1 - time_samples1(1)) / SAMPLING_RATE;
    azimuth_rad1 = deg2rad(azimuth_deg1);
    elevation_rad1 = deg2rad(elevation_deg1);
    R_proj1 = ASSUMED_HEIGHT ./ tan(elevation_rad1);
    x_coords1 = R_proj1 .* cos(azimuth_rad1);
    y_coords1 = R_proj1 .* sin(azimuth_rad1);
    valid_indices1 = isfinite(x_coords1) & isfinite(y_coords1);
    time_sec1 = time_sec1(valid_indices1);
    x_coords1 = x_coords1(valid_indices1);
    y_coords1 = y_coords1(valid_indices1);
    
    % 计算事件1速度
    [velocities1, time_midpoints1, avg_vel1, ~] = calculate_velocity(time_sec1, x_coords1, y_coords1, TIME_WINDOW_VEL);
end

% --- 处理事件2 (左子图数据) ---
fprintf('\n--- 正在处理事件2 (%.2e to %.2e) ---\n', event2_range(1), event2_range(2));
velocities2 = []; time_midpoints2 = []; avg_vel2 = NaN; % 初始化结果
% 筛选事件2数据
logicalIndex2 = ...
    abs(data_raw(:, COL_RCORR)) > 0.6 & ...
    data_raw(:, COL_START_LOC) > event2_range(1) & ...
    data_raw(:, COL_START_LOC) < event2_range(2) & ...
    data_raw(:, COL_ELEVATION) < 80 & ...
    abs(data_raw(:, COL_T123)) < 1;
data2 = data_raw(logicalIndex2, :);
fprintf('事件2筛选后有效数据点数: %d\n', size(data2, 1));

if ~isempty(data2)
    % 预处理事件2数据
    time_samples2 = data2(:, COL_START_LOC);
    azimuth_deg2 = data2(:, COL_AZIMUTH);
    elevation_deg2 = data2(:, COL_ELEVATION);
    time_sec2 = (time_samples2 - time_samples2(1)) / SAMPLING_RATE;
    azimuth_rad2 = deg2rad(azimuth_deg2);
    elevation_rad2 = deg2rad(elevation_deg2);
    R_proj2 = ASSUMED_HEIGHT ./ tan(elevation_rad2);
    x_coords2 = R_proj2 .* cos(azimuth_rad2);
    y_coords2 = R_proj2 .* sin(azimuth_rad2);
    valid_indices2 = isfinite(x_coords2) & isfinite(y_coords2);
    time_sec2 = time_sec2(valid_indices2);
    x_coords2 = x_coords2(valid_indices2);
    y_coords2 = y_coords2(valid_indices2);

    % 计算事件2速度
    [velocities2, time_midpoints2, avg_vel2, ~] = calculate_velocity(time_sec2, x_coords2, y_coords2, TIME_WINDOW_VEL);
end


%% ==================== 3. 结果可视化 ====================
fprintf('\n正在生成对比图...\n');
figure('Name', '双事件速率对比分析', 'Position', [100, 100, 1200, 500]);

% --- 绘制左子图 (事件2) ---
subplot(1, 2, 1);
if ~isempty(velocities2)
    time_plot2 = time_midpoints2 * 1000;
    velo_plot2 = velocities2 / 1e5;
    avg_velo_plot2 = avg_vel2 / 1e5;
    smoothed_velo_plot2 = movmean(velo_plot2, 5);
    
    hold on;
    plot(time_plot2, velo_plot2, 'Color', [0.7 0.7 0.7], 'LineWidth', 1, 'DisplayName', '二维速率');
    plot(time_plot2, smoothed_velo_plot2, 'k', 'LineWidth', 2, 'DisplayName', '五点平滑后');
    line([min(time_plot2), max(time_plot2)], [avg_velo_plot2, avg_velo_plot2], ...
         'Color', 'k', 'LineStyle', '--', 'LineWidth', 1, 'DisplayName', '二维平均速率');
    hold off;
    grid on;
    xlabel('时间/ms');
    ylabel('二维速率/(10^5 m·s^{-1})');
    legend('Location', 'southoutside', 'Orientation', 'horizontal');
    set(gca, 'FontSize', 12);
    xlim([min(time_plot2), max(time_plot2)]);
    ylim([0, 20]);
else
    title(sprintf('事件 2 (%.2e - %.2e)\n无足够数据点', event2_range(1), event2_range(2)));
end

% --- 绘制右子图 (事件1) ---
subplot(1, 2, 2);
if ~isempty(velocities1)
    time_plot1 = time_midpoints1 * 1000;
    velo_plot1 = velocities1 / 1e5;
    avg_velo_plot1 = avg_vel1 / 1e5;
    smoothed_velo_plot1 = movmean(velo_plot1, 5);

    hold on;
    plot(time_plot1, velo_plot1, 'Color', [0.7 0.7 0.7], 'LineWidth', 1, 'DisplayName', '二维速率');
    plot(time_plot1, smoothed_velo_plot1, 'k', 'LineWidth', 2, 'DisplayName', '五点平滑后');
    line([min(time_plot1), max(time_plot1)], [avg_velo_plot1, avg_velo_plot1], ...
         'Color', 'k', 'LineStyle', '--', 'LineWidth', 1, 'DisplayName', '二维平均速率');
    hold off;
    grid on;
    xlabel('时间/ms');
    ylabel('二维速率/(10^5 m·s^{-1})');
    legend('Location', 'southoutside', 'Orientation', 'horizontal');
    set(gca, 'FontSize', 12);
    xlim([min(time_plot1), max(time_plot1)]);
    ylim([0, 10]);
else
    title(sprintf('事件 1 (%.2e - %.2e)\n无足够数据点', event1_range(1), event1_range(2)));
end

fprintf('绘图完成。\n');


%% ==================== 函数定义区 ====================
function [velocities, time_midpoints, avg_vel, std_vel] = calculate_velocity(t, x, y, time_window)
% calculate_velocity: 计算先导在二维平面上的发展速度
% 【修改】: 不再输出高度信息
    velocities = []; time_midpoints = [];
    total_duration = t(end) - t(1);
    num_segments = floor(total_duration / time_window);
    if num_segments == 0
        warning('总时长小于一个时间窗口，无法进行分段速度计算。');
        avg_vel = NaN; std_vel = NaN; return;
    end
    for i = 1:num_segments
        t_start = t(1) + (i-1) * time_window;
        t_end = t_start + time_window;
        indices_in_window = find(t >= t_start & t < t_end);
        if length(indices_in_window) >= 2
            x_segment = x(indices_in_window);
            y_segment = y(indices_in_window);
            t_segment = t(indices_in_window);
            path_length = sum(sqrt(diff(x_segment).^2 + diff(y_segment).^2));
            delta_t = t_segment(end) - t_segment(1);
            if delta_t > 0
                velocities = [velocities; path_length / delta_t];
                time_midpoints = [time_midpoints; t_start + time_window/2];
            end
        end
    end
    if ~isempty(velocities)
        avg_vel = mean(velocities);
        std_vel = std(velocities);
    else
        warning('所有时间段内的数据点均少于2个，无法计算任何速度。');
        avg_vel = NaN; std_vel = NaN;
    end
end