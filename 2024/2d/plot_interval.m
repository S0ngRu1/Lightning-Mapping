%% ==================== 0. 清理与参数设置 ====================
clear; 
clc; 
close all;

% --- 用户可调参数 ---
filename = 'result_yld_3.5e8_4e8_window_512_128_去零飘_滤波_加窗_阈值15_30_80.txt';
start_loc = 3.65e8;
end_loc = 3.8e8;
thea = 3000;
num_subplots = 8;

%% ==================== 1. 数据加载、筛选与梯级识别 ====================
fprintf('--- 正在加载并筛选二维结果 ---\n');
try
    result_table = readtable(filename);
catch
    error('数据文件 "%s" 加载失败，请检查文件名或路径。', filename);
end
logicalIndex = abs(result_table.t123) < 1 & abs(result_table.Rcorr) > 0.65 & result_table.Start_loc >= start_loc & result_table.Start_loc < end_loc;
filtered_table = result_table(logicalIndex, :);
if isempty(filtered_table), error('在指定的时间范围内没有找到满足初始条件的数据点。'); end
filtered_table = sortrows(filtered_table, 'Start_loc');
time_diffs = diff(filtered_table.Start_loc);
gap_indices = find(time_diffs > thea);
step_start_indices = [1; gap_indices + 1];
step_end_indices = [gap_indices; height(filtered_table)];
num_total_steps = numel(step_start_indices);
fprintf('共识别出 %d 个梯级。\n', num_total_steps);


%% ==================== 2. 【全新风格】分阶段、分颜色累积可视化 ====================
fprintf('--- 正在生成 %d 个子图进行分色累积可视化 ---\n', num_subplots);

% 计算每个子图应该“新增”多少个梯级
steps_per_subplot = ceil(num_total_steps / num_subplots);

% 【新增】定义颜色
% 为每个新增阶段定义一种独特的、鲜艳的颜色
new_stage_colors = lines(num_subplots); 
% 为已经画过的“旧”阶段定义一个统一的灰色
previous_stage_color = [0.75 0.75 0.75]; 

% 创建一个白色背景的 figure
figure('Color', 'w');
sgtitle(sprintf('负先导分阶段累积发展图 (%.3e - %.3e)', start_loc, end_loc), ...
        'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k');

% 循环创建8个子图
for i = 1:num_subplots
    subplot(2, 4, i);
    hold on;
    
    % a. 确定当前子图“新增”的梯级范围
    start_step_new = (i-1) * steps_per_subplot + 1;
    end_step_new = min(i * steps_per_subplot, num_total_steps);
    
    if start_step_new > num_total_steps
        axis off; continue;
    end
    
    % b. 绘制“旧”数据 (如果不是第一个子图)
    if i > 1
        % 找到之前所有阶段的结束行号
        end_row_previous = step_end_indices(start_step_new - 1);
        previous_data = filtered_table(1:end_row_previous, :);
        % 用灰色绘制旧数据
        scatter(previous_data.Azimuth, previous_data.Elevation, 3, previous_stage_color, 'filled');
    end
    
    % c. 绘制“新”数据
    % 找到当前新增阶段的起止行号
    start_row_new = step_start_indices(start_step_new);
    end_row_new = step_end_indices(end_step_new);
    new_data = filtered_table(start_row_new:end_row_new, :);
    
    % 用当前阶段的独特颜色绘制新数据
    current_stage_color = new_stage_colors(i, :);
    scatter(new_data.Azimuth, new_data.Elevation, 3, current_stage_color, 'filled');
    
    % d. 设置子图样式 (白色背景，黑色文字)
    xlabel('方位角 (Azimuth / °)', 'Color', 'k');
    ylabel('仰角 (Elevation / °)', 'Color', 'k');
    ylim([0 60])
    xlim([160 210])
    set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'Box', 'on');
    grid on;
    
    hold off;
end

% 【新增】手动创建图例来解释颜色含义
% 我们在一个不可见的新坐标轴上手动创建图例
ax_legend = axes('Position', [0, 0, 0.1, 0.1], 'Visible', 'off');
hold(ax_legend, 'on');
legend_handles = gobjects(num_subplots, 1);
for i = 1:num_subplots
    legend_handles(i) = scatter(ax_legend, NaN, NaN, 100, new_stage_colors(i,:), 'filled');
end
legend(legend_handles, ...
    { '阶段1', '阶段2', '阶段3', '阶段4', '阶段5', '阶段6', '阶段7', '阶段8'}, ...
    'Position', [0.91, 0.3, 0.08, 0.4], 'FontSize', 10);
hold(ax_legend, 'off');


disp('绘图完成!');