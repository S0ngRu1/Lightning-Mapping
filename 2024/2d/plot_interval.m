%% ==================== 0. 清理与参数设置 ====================
clear; 
clc; 
close all;

% --- 假设 all_S_results 和 all_match_results 已经存在于您的工作区 ---
% 如果没有，请先加载包含这两个变量的 .mat 文件
% load('your_data_file.mat');

% --- 用户可调参数 ---
% 分析的时间范围 (单位: 采样点)
start_loc = 3.965e8;
end_loc = 3.98e8;

% 梯级识别的时间阈值 (单位: 采样点)
thea = 3000;

% DBSCAN 聚类参数
epsilon = 200;      % 邻域半径 (米)
min_points = 3;     % 形成核心簇的最小点数

% 最终绘图的子图数量
num_subplots = 8;

%% ==================== 1. 数据筛选与梯级识别 ====================
fprintf('--- 正在筛选 %.3e 至 %.3e 范围内的数据 ---\n', start_loc, end_loc);

% a. 根据初始条件和时间范围筛选数据
conditions = ([all_match_results.yld_start_loc] >= start_loc) & ...
             ([all_match_results.yld_start_loc] < end_loc) & ...
             ([all_match_results.dlta] < 20000) & ...
             ([all_match_results.r_gccs] > 0.1) & ...
             (abs([all_match_results.R3_value]) < 10000);
         
filtered_indices = find(conditions);

if isempty(filtered_indices)
    error('在指定的时间范围内没有找到满足初始条件的数据点。');
end

% 提取对应的3D坐标和时间戳
filtered_S = all_S_results(filtered_indices, :);
filtered_match_result = all_match_results(filtered_indices);
time_samples = [filtered_match_result.yld_start_loc]';

% b. 按时间排序
[time_samples_sorted, sort_order] = sort(time_samples);
S_sorted = filtered_S(sort_order, :);

% c. 识别梯级
fprintf('--- 正在识别所有梯级 ---\n');
time_diffs = diff(time_samples_sorted);
gap_indices = find(time_diffs > thea);

step_start_indices = [1; gap_indices + 1];
step_end_indices = [gap_indices; length(time_samples_sorted)];
num_total_steps = numel(step_start_indices);
fprintf('共识别出 %d 个梯级。\n', num_total_steps);

%% ==================== 2. 对每个梯级进行DBSCAN去噪 ====================
fprintf('--- 正在对每个梯级进行DBSCAN聚类去噪 ---\n');

cleaned_steps_data = cell(num_total_steps, 1); % 使用cell数组存储清洗后的数据

for i = 1:num_total_steps
    % 提取当前梯级的所有点
    step_idx_range = step_start_indices(i):step_end_indices(i);
    
    if numel(step_idx_range) < min_points
        continue; % 点数太少，跳过
    end
    
    points_in_step = S_sorted(step_idx_range, :);
    times_in_step = time_samples_sorted(step_idx_range);
    
    % 执行DBSCAN
    cluster_labels = dbscan(points_in_step, epsilon, min_points);
    
    % 找到最大的簇作为主路径
    main_cluster_id = mode(cluster_labels(cluster_labels > 0));
    
    if isempty(main_cluster_id)
        continue; % 没有找到主簇
    end
    
    main_path_mask = (cluster_labels == main_cluster_id);
    
    % 存储清洗后的主路径数据 [x, y, z, time]
    cleaned_steps_data{i} = [points_in_step(main_path_mask, :), times_in_step(main_path_mask)];
end

% 移除处理后为空的cell
cleaned_steps_data = cleaned_steps_data(~cellfun('isempty', cleaned_steps_data));
num_cleaned_steps = numel(cleaned_steps_data);
fprintf('聚类后剩余 %d 个有效梯级。\n', num_cleaned_steps);

%% ==================== 3. 分组并进行可视化 ====================
fprintf('--- 正在生成 %d 个子图进行可视化 ---\n', num_subplots);

% 计算每个子图应该包含多少个梯级
steps_per_subplot = ceil(num_cleaned_steps / num_subplots);

% 创建一个深色背景的 figure
figure('Color', [0.1 0.1 0.2], 'Position', [50, 50, 1600, 800]);
sgtitle(sprintf('负先导分梯级发展图 (%.3e - %.3e)', start_loc, end_loc), ...
        'FontSize', 18, 'FontWeight', 'bold', 'Color', 'w');

% 循环创建8个子图
for i = 1:num_subplots
    subplot(2, 4, i); % 使用2x4的布局
    hold on;
    
    % a. 确定当前子图要绘制的梯级范围
    start_step_idx = (i-1) * steps_per_subplot + 1;
    end_step_idx = min(i * steps_per_subplot, num_cleaned_steps);
    
    if start_step_idx > num_cleaned_steps
        % 如果没有更多梯级可画，则创建一个空坐标轴
        set(gca, 'Color', [0.1 0.1 0.2], 'XColor', 'w', 'YColor', 'w');
        title(sprintf('梯级 %d - %d', start_step_idx, end_step_idx), 'Color', 'w');
        axis off; % 关闭坐标轴
        continue;
    end
    
    % b. 合并当前子图要绘制的所有点
    subplot_data = vertcat(cleaned_steps_data{start_step_idx:end_step_idx});
    
    x_plot = subplot_data(:, 1);
    y_plot = subplot_data(:, 2);
    time_plot = subplot_data(:, 4);
    
    % c. 归一化时间并绘图 (使用您提供的风格)
    colorValues = (time_plot - time_samples_sorted(1)) / (time_samples_sorted(end) - time_samples_sorted(1));
    
    scatter(x_plot, y_plot, 10, colorValues, 'filled', 'MarkerFaceAlpha', 0.8);
    
    % d. 设置子图样式
    title(sprintf('梯级 %d - %d', start_step_idx, end_step_idx), 'FontSize', 12, 'Color', 'w');
    xlabel('X (m)', 'Color', 'w');
    ylabel('Y (m)', 'Color', 'w');
    set(gca, 'Color', [0.1 0.1 0.2], 'XColor', 'w', 'YColor', 'w', 'GridLineStyle', '--', 'GridAlpha', 0.4, 'Box', 'on');
    grid on;
    axis equal; % 保持XY轴比例一致
    
    hold off;
end

% 为整个 figure 添加一个颜色条
h = colorbar('Position', [0.92 0.1 0.015 0.8]); % 手动指定颜色条位置
colormap('parula');
ylabel(h, '归一化发展时间', 'FontSize', 11, 'Color', 'w');
set(h, 'Color', 'w');
caxis([0, 1]);

disp('绘图完成!');