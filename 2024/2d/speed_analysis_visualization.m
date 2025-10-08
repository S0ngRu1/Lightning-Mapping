%% ==================== 1. 数据准备 ====================
x_range = [-10000, 6000];
y_range = [-10000, 0];
z_range = [0, 10000];
% 筛选条件
conditions = ([all_match_results.dlta] < 20000) & ...
             ([all_match_results.yld_start_loc] > 3.88e8) & ...
             ([all_match_results.yld_start_loc] < 4e8) & ...
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
NUM_SEGMENTS = 400;  
EPSILON = 300;      % 邻域半径设为500米。如果您的通道发展很密集，可以减小此值
MIN_POINTS = 2;     % 至少4个点才能构成一个核心簇
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
[velocities, ~, avg_heights] = calculate_3d_velocity_by_points( ...
    time_sec_sorted, x_coords_sorted, y_coords_sorted, z_coords_sorted, ...
    NUM_SEGMENTS, EPSILON, MIN_POINTS);

%% ==================== 3. 结果可视化 ====================
fprintf('正在生成分析图...\n');
% --- 绘制左子图 (速度 vs. 时间) ---
figure
subplot(1, 2, 1);
if ~isempty(velocities)
    % 准备绘图数据
    time_plot_ms = linspace(0,50,length(velocities));
    velo_plot_1e6 = velocities / 1e7;     % 速度单位: 10^5 m/s
    avg_velo_plot = mean(velo_plot_1e6, 'omitnan');
    smoothed_velo_plot = movmean(velo_plot_1e6, 4, 'omitnan'); % 5点滑动平均
    
    hold on;
    % 绘制原始速度、平滑速度和平均速度
    plot(time_plot_ms, velo_plot_1e6, 'Color', [0.8 0.8 0.8], 'LineWidth', 1, 'DisplayName', '瞬时三维速率');
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
    xlim([0 50]);
    xticks(0:10:50);
    ylim([0 15]);
    yticks(0:3:15);
else
    title('无足够数据进行速度-时间分析');
end

ax_right = subplot(1, 2, 2); % 获取右子图的句柄

% ====================【新增代码】====================
% 获取当前子图的默认位置
pos = get(ax_right, 'Position');
% 修改位置参数：将左边距向右移动一点，同时将宽度缩减
pos(1) = pos(1) + 0.05; % 'left'值增加，整体右移
pos(3) = pos(3) * 0.75; % 'width'值缩小为原来的75%
% 应用新的位置设置
set(ax_right, 'Position', pos);
% ===================================================

if ~isempty(velocities)
    % 准备绘图数据
    height_plot_km = avg_heights; % 高度是米，无需转换
    
    % 使用散点图，并用颜色表示时间
    scatter(velo_plot_1e6, height_plot_km, 30, time_plot_ms, "filled");
    
    grid on;
    xlabel('三维速率 (10^5 m/s)');
    ylabel('高度 (m)');
    title('速率与放电高度的关系');
    xlim([0 10]);
    xticks(0:2:10);
    % 添加颜色条并标注
    h_bar = colorbar;
    ylabel(h_bar, '时间 (ms)'); % YLabel的单位应与颜色数据time_plot_ms一致
    colormap('jet'); % 使用jet颜色图
    
    set(gca, 'FontSize', 12);
else
    title('无足够数据进行速度-高度分析');
end

%% ==================== 函数定义区 ====================
function [velocities, time_midpoints, avg_heights] = calculate_3d_velocity_by_points(t, x, y, z, num_segments, epsilon, min_points)
% calculate_3d_velocity_by_points: 
%   将总数据点分成 num_segments 份，对每一份数据进行聚类和速度计算。
%
%   num_segments: 您希望将总数据点分成多少份来进行计算 (例如 500)。

    velocities = [];
    time_midpoints = [];
    avg_heights = [];
    
    total_points = length(t);
    if total_points < min_points
        warning('总数据点数太少，无法进行计算。');
        return;
    end
    
    % 计算每一份(段)应该包含多少个点
    points_per_segment = floor(total_points / num_segments);
    
    for i = 1:num_segments
        % 计算当前段的起始和结束索引
        start_idx = (i-1) * points_per_segment + 1;
        end_idx = i * points_per_segment;
        
        % 确保最后一个段包含所有剩余的点
        if i == num_segments
            end_idx = total_points;
        end
        
        indices_in_window = start_idx:end_idx;
        
        if length(indices_in_window) >= min_points
            % --- 后续的聚类和路径计算逻辑与之前完全相同 ---
            
            % 1. 在点窗内进行DBSCAN聚类
            points_in_window = [x(indices_in_window), y(indices_in_window), z(indices_in_window)];
            idx = dbscan(points_in_window, epsilon, min_points);
            main_cluster_id = mode(idx(idx > 0));
            
            if isempty(main_cluster_id)
                continue;
            end
            
            % 2. 筛选出属于主路径的点
            main_path_mask = (idx == main_cluster_id);
            if sum(main_path_mask) < 2
                continue;
            end
            
            x_main_path = points_in_window(main_path_mask, 1);
            y_main_path = points_in_window(main_path_mask, 2);
            z_main_path = points_in_window(main_path_mask, 3);
            t_main_path = t(indices_in_window(main_path_mask));
            
            [t_main_path_sorted, sort_order] = sort(t_main_path);
            x_main_path_sorted = x_main_path(sort_order);
            y_main_path_sorted = y_main_path(sort_order);
            z_main_path_sorted = z_main_path(sort_order);

            % 3. 仅对主路径计算速度
            path_length = sum(sqrt(diff(x_main_path_sorted).^2 + diff(y_main_path_sorted).^2 ));         
            delta_t = t_main_path_sorted(end) - t_main_path_sorted(1);
            
            if delta_t > 0
                velocities = [velocities; path_length / delta_t];
                % 时间中点使用当前段内主路径的平均时间，更准确
                time_midpoints = [time_midpoints; mean(t_main_path_sorted)];
                avg_heights = [avg_heights; mean(z_main_path_sorted)];
            end
        end
    end
    
    if isempty(velocities)
        warning('未能计算出任何有效的速度值。');
    end
end