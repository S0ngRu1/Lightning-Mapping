%% ==================== 1. 数据准备 ====================
%% 新版数据结构
all_match_results = readtable('results\3d_win512_cost_cal_yld_chj_dtoa3.6e8_4.0e8.csv');
% 筛选条件
conditions = ([all_match_results.dlta] < 20000) & ...
             ([all_match_results.yld_start_loc] > 3.82e8) & ...
             ([all_match_results.yld_start_loc] < 3.9e8) & ...
             ([all_match_results.x] > -10000) & ...
             ([all_match_results.x] < 6000) & ...
             ([all_match_results.y] > -10000) & ...
             ([all_match_results.y] < 0) & ...
             ([all_match_results.z] > 0) & ...
             ([all_match_results.z] < 10000) & ...
             ([all_match_results.r_gccs] > 0.1) & ...
             (abs([all_match_results.R3_value]) < 10000);
filtered_match_indices = find(conditions);
filtered_match_result = all_match_results(filtered_match_indices, :);

% 定义时间区间和对应的EPSILON
time_intervals = [
    3.655e8,  3.66e8;   % 区间1
    3.66e8,   3.67e8;   % 区间2
    3.68e8,   3.7e8;    % 区间3
    3.7e8,    3.72e8;   % 区间4
    3.82e8,   3.87e8    % 区间5
];
EPSILON = [45, 100, 120, 100, 45];  % 对应5个区间的EPSILON

% 为每个数据点分配对应的EPSILON（根据其yld_start_loc所在区间）
epsilon_values = zeros(size(filtered_match_result.yld_start_loc));  % 存储每个点的EPSILON
for k = 1:length(filtered_match_result.yld_start_loc)
    yld = filtered_match_result.yld_start_loc(k);
    % 判断当前点属于哪个时间区间
    for m = 1:size(time_intervals, 1)
        if yld >= time_intervals(m, 1) && yld < time_intervals(m, 2)
            epsilon_values(k) = EPSILON(m);
            break;  % 找到对应区间后退出内层循环
        end
    end
end

% 过滤不在任何时间区间内的点（如果有的话）
invalid_idx = epsilon_values == 0;
if any(invalid_idx)
    warning('有%d个点不在指定的时间区间内，已自动过滤。', sum(invalid_idx));
    filtered_match_result = filtered_match_result(~invalid_idx, :);
    epsilon_values = epsilon_values(~invalid_idx);
end

% 从数据结构中提取时间和坐标（过滤后的数据）
time_samples = [filtered_match_result.yld_start_loc]';
x_coords = [filtered_match_result.x];
y_coords = [filtered_match_result.y]; 
z_coords = [filtered_match_result.z]; 

% --- 用户可调参数 ---
SAMPLING_RATE = 200e6;       % 数据采集卡采样率 (Hz), 200 MS/s
NUM_SEGMENTS = 40;  
MIN_POINTS = 4;     % 至少4个点才能构成一个核心簇

%% ==================== 2. 数据预处理和速度计算 ====================

% 确保数据点数大于1
if numel(time_samples) < 2
    error('有效数据点不足2个，无法计算速度。');
end

% 将采样点时间转换为以秒为单位的相对时间
time_sec = (time_samples - 0) / SAMPLING_RATE;

% 调用函数计算三维速度（传入每个点的EPSILON值）
fprintf('正在计算三维发展速度...\n');
[velocities, ~, avg_heights] = calculate_3d_velocity_by_points( ...
    time_sec, x_coords, y_coords, z_coords, ...
    NUM_SEGMENTS, epsilon_values, MIN_POINTS);

valid_indices = velocities <= 2e6;
velocities = velocities(valid_indices);
avg_heights = avg_heights(valid_indices);

%% ==================== 3. 结果可视化 ====================
fprintf('正在生成分析图...\n');
% --- 绘制左子图 (速度 vs. 时间) ---
figure
subplot(1, 2, 1);
if ~isempty(velocities)
    % 准备绘图数据
    time_plot_ms = linspace(0,70,length(velocities));
    velo_plot_1e5 = velocities / 1e5;     % 速度单位: 10^5 m/s
    avg_velo_plot = mean(velo_plot_1e5, 'omitnan');
    smoothed_velo_plot = movmean(velo_plot_1e5, 2, 'omitnan'); % 滑动平均
    
    hold on;
    % 绘制原始速度、平滑速度和平均速度
    plot(time_plot_ms, velo_plot_1e5, 'Color', [0.8 0.8 0.8], 'LineWidth', 1, 'DisplayName', '瞬时三维速率');
    plot(time_plot_ms, smoothed_velo_plot, 'b-', 'LineWidth', 1, 'DisplayName', '平滑后速率');
    line([min(time_plot_ms), max(time_plot_ms)], [avg_velo_plot, avg_velo_plot], ...
         'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.5, 'DisplayName', '平均速率');
    hold off;
    
    grid on;
    xlabel('时间 (ms)');
    ylabel('三维速率 (10^5 m/s)');
    title('闪电发展速率随时间的变化');
    legend('show', 'Location', 'northwest');
    set(gca, 'FontSize', 12);
    xlim([0 40]);
    xticks(0:5:40);
    % ylim([0 15]);
    % yticks(0:3:15);
else
    title('无足够数据进行速度-时间分析');
end

ax_right = subplot(1, 2, 2); % 获取右子图的句柄

% 调整右子图位置
pos = get(ax_right, 'Position');
pos(1) = pos(1) + 0.05; % 右移
pos(3) = pos(3) * 0.75; % 缩减宽度
set(ax_right, 'Position', pos);

if ~isempty(velocities)
    % 准备绘图数据
    height_plot_km = avg_heights; % 高度单位：米
    
    % 散点图，颜色表示时间
    scatter(velo_plot_1e5, height_plot_km, 30, time_plot_ms, "filled");
    
    grid on;
    xlabel('三维速率 (10^5 m/s)');
    ylabel('高度 (m)');
    title('速率与放电高度的关系');
    % 添加颜色条
    h_bar = colorbar;
    ylabel(h_bar, '时间 (ms)');
    colormap('jet');
    
    set(gca, 'FontSize', 12);
else
    title('无足够数据进行速度-高度分析');
end

%% ==================== 函数定义区 ====================
function [velocities, time_midpoints, avg_heights] = calculate_3d_velocity_by_points(t, x, y, z, num_segments, epsilon_values, min_points)
% calculate_3d_velocity_by_points: 
%   将总数据点分成 num_segments 份，对每一份数据根据时间区间对应的EPSILON进行聚类和速度计算。

    velocities = [];
    time_midpoints = [];
    avg_heights = [];
    
    total_points = length(t);
    if total_points < min_points
        warning('总数据点数太少，无法进行计算。');
        return;
    end
    
    % 计算每段包含的点数
    points_per_segment = floor(total_points / num_segments);
    
    for i = 1:num_segments
        % 计算当前段的起始和结束索引
        start_idx = (i-1) * points_per_segment + 1;
        end_idx = i * points_per_segment;
        
        % 最后一段包含剩余所有点
        if i == num_segments
            end_idx = total_points;
        end
        
        indices_in_window = start_idx:end_idx;
        
        if length(indices_in_window) >= min_points
            % 获取当前段内所有点的EPSILON值，取众数作为该段的聚类参数
            segment_epsilons = epsilon_values(indices_in_window);
            segment_epsilon = mode(segment_epsilons);  % 以段内多数点的EPSILON为准
            
            % 1. 在当前段内进行DBSCAN聚类（使用该段的EPSILON）
            points_in_window = [x(indices_in_window), y(indices_in_window), z(indices_in_window)];
            idx = dbscan(points_in_window, segment_epsilon, min_points);
            main_cluster_id = mode(idx(idx > 0));  % 主聚类ID
            
            if isempty(main_cluster_id)
                continue;
            end
            
            % 2. 筛选主路径点
            main_path_mask = (idx == main_cluster_id);
            if sum(main_path_mask) < 2  % 至少需要2个点计算速度
                continue;
            end
            
            x_main_path = points_in_window(main_path_mask, 1);
            y_main_path = points_in_window(main_path_mask, 2);
            z_main_path = points_in_window(main_path_mask, 3);
            t_main_path = t(indices_in_window(main_path_mask));
            
            % 按时间排序
            [t_main_path_sorted, sort_order] = sort(t_main_path);
            x_main_path_sorted = x_main_path(sort_order);
            y_main_path_sorted = y_main_path(sort_order);
            z_main_path_sorted = z_main_path(sort_order);

            % 3. 计算速度
            path_length = sum(sqrt(diff(x_main_path_sorted).^2 + diff(y_main_path_sorted).^2 ));         
            delta_t = t_main_path_sorted(end) - t_main_path_sorted(1);
            
            if delta_t > 0
                velocities = [velocities; path_length / delta_t];
                time_midpoints = [time_midpoints; mean(t_main_path_sorted)];
                avg_heights = [avg_heights; mean(z_main_path_sorted)];
            end
        end
    end
    
    if isempty(velocities)
        warning('未能计算出任何有效的速度值。');
    end
end