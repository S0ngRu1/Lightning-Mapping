%% ==================== 1. 数据准备与参数定义 ====================
all_match_results = readtable('results\3d_win512_cost_cal_yld_chj_dtoa3.6e8_4.0e8.csv');
% --- 定义您要分析的六个时间段 ---
time_intervals = [
    3.65e8,  3.66e8;
    3.66e8, 3.67e8;
    3.67e8,  3.7e8;
    3.7e8, 3.72e8;
    3.82e8,  3.9e8
    ];

% --- 用户可调参数 ---
SAMPLING_RATE = 200e6;
NUM_SEGMENTS_PER_INTERVAL = 1; % 在每个时间段内，再细分成的计算单元数
EPSILON = 80;
MIN_POINTS = 2;

fprintf('开始计算每个时间段的平均速度...\n');
fprintf('----------------------------------------\n');

%% ==================== 2. 循环计算每个时间段的速度 ====================

for i = 1:size(time_intervals, 1)

    current_start_loc = time_intervals(i, 1);
    current_end_loc = time_intervals(i, 2);

    fprintf('正在处理时间段: %.3e to %.3e ...\n', current_start_loc, current_end_loc);

    % 筛选条件
    conditions = ([all_match_results.dlta] < 6000) & ...
        ([all_match_results.yld_start_loc] > current_start_loc) & ...
        ([all_match_results.yld_start_loc] < current_end_loc) & ...
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
    % 从数据结构中提取时间和坐标
    time_samples = [filtered_match_result.yld_start_loc]';
    x_coords= [filtered_match_result.x];
    y_coords = [filtered_match_result.y];
    z_coords = [filtered_match_result.z];
    %% ==================== 2. 数据预处理和速度计算 ====================

    % 确保数据点数大于1
    if numel(time_samples) < 2
        error('有效数据点不足2个，无法计算速度。');
    end

    % 将采样点时间转换为以秒为单位的相对时间
    time_sec = (time_samples - 0) / SAMPLING_RATE;

    % 调用函数计算速度
    [velocities, ~, ~] = calculate_3d_velocity_by_points( ...
        time_sec, x_coords, y_coords, z_coords, ...
        NUM_SEGMENTS_PER_INTERVAL, EPSILON, MIN_POINTS);

    % --- c. 计算平均速度并输出结果 ---
    if ~isempty(velocities)
        avg_speed = mean(velocities, 'omitnan');
        if avg_speed > 1e6
            fprintf('  -> 平均速度为: %.2f m/s\n\n', avg_speed/1e6);
        else
            fprintf('  -> 平均速度为: %.2f m/s\n\n', avg_speed/1e5);
        end
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
        path_length = sum(sqrt(diff(x_main_path_sorted).^2 + diff(y_main_path_sorted).^2 ));
        delta_t = t_main_path_sorted(end) - t_main_path_sorted(1);

        if delta_t > 0
            velocities = [velocities; path_length / delta_t];
            time_midpoints = [time_midpoints; mean(t_main_path_sorted)];
            avg_heights = [avg_heights; mean(z_main_path_sorted)];
        end
    end
end
end