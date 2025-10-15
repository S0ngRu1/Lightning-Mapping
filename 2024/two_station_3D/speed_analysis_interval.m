%% ==================== 1. 数据准备与参数定义 ====================
% --- 假设 all_S_results 和 all_match_results 已经存在于您的工作区 ---
% ...

% --- 定义您要分析的六个时间段 ---
time_intervals = [
    3.8e8,  3.84e8;
    3.84e8, 3.88e8;
    3.989e8, 4e8;
    3.965e8, 3.98e8;
    3.94e8, 3.965e8;
    3.9e8,  3.94e8
];

% --- 用户可调参数 ---
x_range = [-10000, 6000];
y_range = [-10000, 0];
z_range = [0, 10000];
SAMPLING_RATE = 200e6;
NUM_SEGMENTS_PER_INTERVAL = 100; % 在每个时间段内，再细分成的计算单元数
EPSILON = 300;
MIN_POINTS = 4;

fprintf('开始计算每个时间段的平均速度...\n');
fprintf('----------------------------------------\n');

%% ==================== 2. 循环计算每个时间段的速度 ====================

for i = 1:size(time_intervals, 1)
    
    current_start_loc = time_intervals(i, 1);
    current_end_loc = time_intervals(i, 2);
    
    fprintf('正在处理时间段: %.3e to %.3e ...\n', current_start_loc, current_end_loc);
    
    % --- a. 针对当前时间段筛选数据 ---
    conditions = ([all_match_results.dlta] < 20000) & ...
                 ([all_match_results.yld_start_loc] >= current_start_loc) & ...
                 ([all_match_results.yld_start_loc] < current_end_loc) & ...
                 ([all_match_results.r_gccs] > 0.1) & ...
                 (abs([all_match_results.R3_value]) < 10000);
             
    filtered_match_indices = find(conditions);
    
    if isempty(filtered_match_indices)
        fprintf('  -> 该时间段内无满足初始条件的数据点。\n\n');
        continue;
    end
    
    filtered_S_temp = all_S_results(filtered_match_indices, :);
    filtered_match_result_temp = all_match_results(filtered_match_indices);
    
    range_condition_s = filtered_S_temp(:,1) >= x_range(1) & filtered_S_temp(:,1) <= x_range(2) & ...
        filtered_S_temp(:,2) >= y_range(1) & filtered_S_temp(:,2) <= y_range(2) & ...
        filtered_S_temp(:,3) >= z_range(1) & filtered_S_temp(:,3) <= z_range(2);
        
    filtered_S = filtered_S_temp(range_condition_s, :);
    filtered_match_result = filtered_match_result_temp(range_condition_s);

    time_samples = [filtered_match_result.yld_start_loc]';
    x_coords = filtered_S(:, 1);
    y_coords = filtered_S(:, 2);
    z_coords = filtered_S(:, 3);
    
    % --- b. 数据预处理和速度计算 ---
    if numel(time_samples) < MIN_POINTS
        fprintf('  -> 该时间段内有效数据点不足 (%d个)，无法计算速度。\n\n', numel(time_samples));
        continue;
    end
    
    time_sec = time_samples / SAMPLING_RATE;
    [time_sec_sorted, sort_indices] = sort(time_sec);
    x_coords_sorted = x_coords(sort_indices);
    y_coords_sorted = y_coords(sort_indices);
    z_coords_sorted = z_coords(sort_indices);
    
    % 调用函数计算速度
    [velocities, ~, ~] = calculate_3d_velocity_by_points( ...
        time_sec_sorted, x_coords_sorted, y_coords_sorted, z_coords_sorted, ...
        NUM_SEGMENTS_PER_INTERVAL, EPSILON, MIN_POINTS);
    
    % --- c. 计算平均速度并输出结果 ---
    if ~isempty(velocities)
        avg_speed = mean(velocities, 'omitnan');
        fprintf('  -> 平均速度为: %.2f m/s\n\n', avg_speed/1e7);
    else
        fprintf('  -> 未能计算出有效速度值。\n\n');
    end
    
end

fprintf('----------------------------------------\n');
fprintf('所有时间段计算完成。\n');


%% ==================== 函数定义区 ====================
function [velocities, time_midpoints, avg_heights] = calculate_3d_velocity_by_points(t, x, y, z, num_segments, epsilon, min_points)
    velocities = [];
    time_midpoints = [];
    avg_heights = [];
    
    total_points = length(t);
    if total_points < min_points
        return;
    end
    
    points_per_segment = floor(total_points / num_segments);
    if points_per_segment == 0
        points_per_segment = 1; % 确保至少有一个点
        num_segments = total_points;
    end
    
    for i = 1:num_segments
        start_idx = (i-1) * points_per_segment + 1;
        end_idx = i * points_per_segment;
        if i == num_segments
            end_idx = total_points;
        end
        indices_in_window = start_idx:end_idx;
        
        if length(indices_in_window) >= min_points
            points_in_window = [x(indices_in_window), y(indices_in_window), z(indices_in_window)];
            idx = dbscan(points_in_window, epsilon, min_points);
            main_cluster_id = mode(idx(idx > 0));
            
            if isempty(main_cluster_id)
                continue;
            end
            
            main_path_mask = (idx == main_cluster_id);
            if sum(main_path_mask) < 2
                continue;
            end
            
            t_main_path = t(indices_in_window(main_path_mask));
            x_main_path = points_in_window(main_path_mask, 1);
            y_main_path = points_in_window(main_path_mask, 2);
            z_main_path = points_in_window(main_path_mask, 3);
            
            [t_main_path_sorted, sort_order] = sort(t_main_path);
            x_main_path_sorted = x_main_path(sort_order);
            y_main_path_sorted = y_main_path(sort_order);
            z_main_path_sorted = z_main_path(sort_order);

            % 【修正】确保计算的是三维路径长度
            path_length = sum(sqrt(diff(x_main_path_sorted).^2 + diff(y_main_path_sorted).^2 + diff(z_main_path_sorted).^2));         
            delta_t = t_main_path_sorted(end) - t_main_path_sorted(1);
            
            if delta_t > 0
                velocities = [velocities; path_length / delta_t];
                time_midpoints = [time_midpoints; mean(t_main_path_sorted)];
                avg_heights = [avg_heights; mean(z_main_path_sorted)];
            end
        end
    end
end