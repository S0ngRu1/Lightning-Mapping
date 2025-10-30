%% ==================== 0. 清理与参数设置 ====================
clear; 
clc; 
close all;
filename = 'results\20240822165932_result_yld_3.6e8_5.6e8_window_1024_256_阈值4倍标准差_去零飘_30_80_hann.txt';
start_loc = 3.711e8;
end_loc = 3.712e8;
thea = 3000;
point_size = 20;  
num_subplots = 20;  
rows = 4;  % 子图行数
cols = 5;  % 子图列数
%% ==================== 1. 数据加载、筛选与梯级识别 ====================
fprintf('--- 正在加载并筛选二维结果 ---\n');
try
    result_table = readtable(filename);
catch
    error('数据文件 "%s" 加载失败，请检查文件名或路径。', filename);
end

% 数据筛选条件
logicalIndex = abs(result_table.t123) < 1 & ...
               abs(result_table.Rcorr) > 0.6 & ...
               result_table.Start_loc >= start_loc & ...
               result_table.Start_loc < end_loc & ...
               result_table.Elevation < 80 ;

filtered_table = result_table(logicalIndex, :);

if isempty(filtered_table)
    error('在指定的时间范围内没有找到满足初始条件的数据点。');
end

% 按时间排序并识别梯级
filtered_table = sortrows(filtered_table, 'Start_loc');
time_diffs = diff(filtered_table.Start_loc);
gap_indices = find(time_diffs > thea);

if length(gap_indices) < 2
    error('找到的间隙数量不足，无法确定有效梯级。');
end

% 确定目标梯级（最大间隙对应的梯级）
gap_duration = diff(gap_indices);
[~, index] = max(gap_duration);
step_start_indice = gap_indices(index);    % 梯级起始索引（间隙位置）
step_end_indice = gap_indices(index + 1);  % 梯级结束索引（间隙位置）

% 提取当前梯级的所有数据点
current_step_data = filtered_table(step_start_indice+1 : step_end_indice, :);
num_points = height(current_step_data);

if num_points == 0
    error('目标梯级内没有数据点，请检查参数设置。');
end

%% ==================== 2. 梯级数据分份与累积可视化（灰色累积+彩色新增） ====================
fprintf('--- 正在生成 %d 个子图进行分色累积可视化 ---\n', num_subplots);

% 定义颜色：新增部分用独特颜色，累积部分用灰色（保持逻辑）
new_stage_colors = lines(num_subplots);  % 20种区分度高的颜色
previous_stage_color = [0.75 0.75 0.75];  % 累积部分的灰色

% 计算每份的点数（最后几份可能略多，保持原有分配逻辑）
points_per_part = floor(num_points / num_subplots);
remainder = num_points - points_per_part * num_subplots;

% 计算每份的起止索引
part_indices = cell(num_subplots, 1);
current_start = 1;
for p = 1:num_subplots
    current_end = current_start + points_per_part - 1;
    if p <= remainder  % 余数分配到前面的份数
        current_end = current_end + 1;
    end
    part_indices{p} = current_start:current_end;
    current_start = current_end + 1;
end

% 创建图形（增大窗口尺寸以容纳20个子图）
figure('Color', 'w', 'Position', [100 100 1400 900]);  % 宽×高适当增大
sgtitle(sprintf('单个梯级分阶段累积图 (%d - %d)', ...
                current_step_data.Start_loc(1), current_step_data.Start_loc(end)), ...
        'FontSize', 16, 'FontWeight', 'bold', 'Color', 'k');

% 循环绘制子图：前i-1份（灰色）+第i份（彩色）
for i = 1:num_subplots
    subplot(rows, cols, i);  % 改为行列布局（4行5列）
    hold on; grid on;
    
    % 1. 绘制前i-1份的累积数据（灰色）
    if i > 1  
        prev_indices = [];
        for j = 1:i-1
            prev_indices = [prev_indices; part_indices{j}'];  % 累积前i-1份
        end
        prev_data = current_step_data(prev_indices, :);
        scatter(prev_data.Azimuth, prev_data.Elevation, point_size, previous_stage_color, 'filled');
    end
    
    % 2. 绘制第i份的新增数据（彩色）
    current_data = current_step_data(part_indices{i}, :);
    scatter(current_data.Azimuth, current_data.Elevation, point_size, new_stage_colors(i,:), 'filled');
    % 图形样式优化（避免标签拥挤）
    % 仅第一列子图显示y轴标签（避免重复）
    if mod(i-1, cols) == 0  % 第一列：i=1,6,11,16
        ylabel('仰角 (°)', 'Color', 'k', 'FontSize', 7);
    end
    % 仅最后一行子图显示x轴标签（避免重复）
    if i > rows*(cols-1)  % 最后一行：i=16~20
        xlabel('方位角 (°)', 'Color', 'k', 'FontSize', 7);
    end

    % 统一坐标轴范围和样式
    % xlim([155 175]);
    % ylim([65 70]);
    set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'Box', 'on', ...
        'FontSize', 6, 'TickDir', 'in');  % 减小字体，刻度向内更整洁
    axis equal;  % 保持坐标比例一致
    hold off;
end

disp('绘图完成!');

fprintf(' %.2f us , 每个图的分辨率为%.2f us \n', (current_step_data.Start_loc(end)-current_step_data.Start_loc(1))*5/1e3,(current_step_data.Start_loc(end)-current_step_data.Start_loc(1))*5/1e3/20);