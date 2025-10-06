%% ==================== 1. 数据准备 ====================
x_range = [-10000, 6000];
y_range = [-10000, 0];
z_range = [0, 10000];
% 筛选条件
conditions = ([all_match_results.dlta] < 20000) & ...
             ([all_match_results.yld_start_loc] > 3.83e8) & ...
             ([all_match_results.yld_start_loc] < 3.89e8) & ...
             ([all_match_results.r_gccs] > 0.1) & ...
             (abs([all_match_results.R3_value]) < 10000);
% 获取满足条件的索引
filtered_match_indices = find(conditions);
filtered_S_temp = all_S_results(filtered_match_indices, :);
filtered_match_result_temp = all_match_results(filtered_match_indices); 

% 范围筛选
range_condition_s = filtered_S_temp(:,1) >= x_range(1) & filtered_S_temp(:,1) <= x_range(2) & ...
    filtered_S_temp(:,2) >= y_range(1) & filtered_S_temp(:,2) <= y_range(2) & ...
    filtered_S_temp(:,3) >= z_range(1) & filtered_S_temp(:,3) <= z_range(2);

filtered_S = filtered_S_temp(range_condition_s, :);
filtered_match_result = filtered_match_result_temp(range_condition_s);

% 从数据结构中提取时间和坐标
time_samples = [filtered_match_result.yld_start_loc]';
x_coords = filtered_S(:, 1);
y_coords = filtered_S(:, 2);
z_coords = filtered_S(:, 3);

% --- 用户可调参数 ---
SAMPLING_RATE = 200e6;       % 数据采集卡采样率 (Hz), 200 MS/s
TIME_WINDOW_VEL = 100e-6;    % 速度计算的时间窗口 (秒), 例如 100 微秒
EPSILON = 500;      % 邻域半径设为500米。如果您的通道发展很密集，可以减小此值
MIN_POINTS = 4;     % 至少4个点才能构成一个核心簇
%% ==================== 2. 数据预处理和速度计算 ====================

% 确保数据点数大于1
if numel(time_samples) < 2
    error('有效数据点不足2个，无法计算速度。');
end

% 将采样点时间转换为以秒为单位的相对时间
time_sec = (time_samples - 0) / SAMPLING_RATE;

% 将所有数据按时间顺序排序
[time_sec_sorted, sort_indices] = sort(time_sec);
x_coords_sorted = x_coords(sort_indices);
y_coords_sorted = y_coords(sort_indices);
z_coords_sorted = z_coords(sort_indices);

% 调用函数计算三维速度
fprintf('正在计算三维发展速度...\n');
[velocities, time_midpoints, avg_heights] = calculate_3d_velocity_clustered( ...
    time_sec_sorted, x_coords_sorted, y_coords_sorted, z_coords_sorted, ...
    TIME_WINDOW_VEL, EPSILON, MIN_POINTS);

%% ==================== 3. 结果可视化 ====================
fprintf('正在生成分析图...\n');
% --- 绘制左子图 (速度 vs. 时间) ---
subplot(1, 2, 1);
if ~isempty(velocities)
    % 准备绘图数据
    time_plot_ms = time_midpoints * 1000; % 时间转换为毫秒
    velo_plot_1e6 = velocities / 1e6;     % 速度单位: 10^5 m/s
    avg_velo_plot = mean(velo_plot_1e6, 'omitnan');
    smoothed_velo_plot = movmean(velo_plot_1e6, 5, 'omitnan'); % 5点滑动平均
    
    hold on;
    % 绘制原始速度、平滑速度和平均速度
    plot(time_plot_ms, velo_plot_1e6, 'Color', [0.8 0.8 0.8], 'LineWidth', 1, 'DisplayName', '瞬时三维速率');
    plot(time_plot_ms, smoothed_velo_plot, 'black-', 'LineWidth', 1, 'DisplayName', '平滑后速率');
    line([min(time_plot_ms), max(time_plot_ms)], [avg_velo_plot, avg_velo_plot], ...
         'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.5, 'DisplayName', '平均速率');
    hold off;
    
    grid on;
    xlabel('时间 (ms)');
    ylabel('三维速率 (10^5 m/s)');
    title('闪电发展速率随时间的变化');
    legend('show', 'Location', 'northwest');
    set(gca, 'FontSize', 12);
    xlim([min(time_plot_ms), max(time_plot_ms)]);
    xticks(0:40:360);
else
    title('无足够数据进行速度-时间分析');
end

% --- 绘制右子图 (速度 vs. 高度) ---
subplot(1, 2, 2);
if ~isempty(velocities)
    % 准备绘图数据
    height_plot_km = avg_heights / 1000; % 高度转换为公里
    
    % 使用散点图，并用颜色表示时间
    scatter(velo_plot_1e6, height_plot_km, 30, time_midpoints, 'filled');
    
    grid on;
    xlabel('三维速率 (10^5 m/s)');
    ylabel('高度 (km)');
    title('速率与放电高度的关系');
    
    % 添加颜色条并标注
    h_bar = colorbar;
    ylabel(h_bar, '时间 (s)');
    colormap('jet'); % 使用jet颜色图
    
    set(gca, 'FontSize', 12);
else
    title('无足够数据进行速度-高度分析');
end

fprintf('绘图完成。\n');

%% ==================== 函数定义区 ====================
function [velocities, time_midpoints, avg_heights] = calculate_3d_velocity_clustered(t, x, y, z, time_window, epsilon, min_points)
% calculate_3d_velocity_clustered: 
% 在滑动时间窗内，首先使用DBSCAN聚类识别主路径，然后计算三维发展速度和平均高度。
%
%   epsilon: DBSCAN的邻域搜索半径（米）。决定了点与点之间多近才算一个簇。
%   min_points: 形成一个簇所需的最少点数。

    velocities = [];
    time_midpoints = [];
    avg_heights = [];
    
    total_duration = t(end) - t(1);
    num_segments = floor(total_duration / time_window);
    
    if num_segments == 0
        warning('总时长小于一个时间窗口，无法进行分段速度计算。');
        return;
    end
    
    for i = 1:num_segments
        t_start = t(1) + (i-1) * time_window;
        t_end = t_start + time_window;
        
        indices_in_window = find(t >= t_start & t < t_end);
        
        % 至少需要 min_points 个点才能进行有意义的聚类和速度计算
        if length(indices_in_window) >= min_points
            
            % 1. 【新增】在时间窗内进行DBSCAN聚类
            points_in_window = [x(indices_in_window), y(indices_in_window), z(indices_in_window)];
            
            % 调用DBSCAN，idx为每个点的簇标签（0代表噪声）
            idx = dbscan(points_in_window, epsilon, min_points);
            
            % 找到最大的那个簇（点数最多的簇）
            % mode函数会返回出现次数最多的元素
            main_cluster_id = mode(idx(idx > 0)); % 只在非噪声点中寻找
            
            % 如果没有找到任何簇，则跳过这个时间窗
            if isempty(main_cluster_id)
                continue;
            end
            
            % 2. 【新增】筛选出属于主路径的点
            main_path_mask = (idx == main_cluster_id);
            
            % 如果主路径上的点数少于2个，也无法计算速度
            if sum(main_path_mask) < 2
                continue;
            end
            
            x_main_path = points_in_window(main_path_mask, 1);
            y_main_path = points_in_window(main_path_mask, 2);
            z_main_path = points_in_window(main_path_mask, 3);
            t_main_path = t(indices_in_window(main_path_mask));
            
            % 【重要】对筛选出的主路径点，必须重新按时间排序
            [t_main_path_sorted, sort_order] = sort(t_main_path);
            x_main_path_sorted = x_main_path(sort_order);
            y_main_path_sorted = y_main_path(sort_order);
            z_main_path_sorted = z_main_path(sort_order);

            % 3. 【修改】仅对主路径计算速度
%             path_length = sum(sqrt(diff(x_main_path_sorted).^2 + diff(y_main_path_sorted).^2 + diff(z_main_path_sorted).^2));
            path_length = sum(sqrt(diff(x_main_path_sorted).^2 + diff(y_main_path_sorted).^2 ));            
delta_t = t_main_path_sorted(end) - t_main_path_sorted(1);
            
            if delta_t > 0
                velocities = [velocities; path_length / delta_t];
                time_midpoints = [time_midpoints; t_start + time_window / 2];
                avg_heights = [avg_heights; mean(z_main_path_sorted)];
            end
        end
    end
    
    if isempty(velocities)
        warning('所有时间段内的主路径数据点均少于2个，无法计算任何速度。');
    end
end