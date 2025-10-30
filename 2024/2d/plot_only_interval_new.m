% 绘制单个梯级的代码，并且先绘制完整数据（灰色作为背景）
%% ==================== 0. 清理与参数设置 ====================
clear; 
clc; 
close all;
filename = 'results\20240822165932_result_yld_3.6e8_5.6e8_window_1024_256_阈值4倍标准差_去零飘_30_80_hann.txt';
step_start_indice = 385909739-1;    % 梯级起始索引（间隙位置）
step_end_indice = 385966721+1;  % 梯级结束索引（间隙位置）
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
               abs(result_table.Rcorr) > 0.4 & ...
               result_table.Start_loc >= step_start_indice & ...
               result_table.Start_loc < step_end_indice & ...
               result_table.Elevation < 80 ;

current_step_data = result_table(logicalIndex, :);

% 提取当前梯级的所有数据点
num_points = height(current_step_data);

if num_points == 0
    error('目标梯级内没有数据点，请检查参数设置。');
end

%% ==================== 2. 梯级数据分份与可视化（完整灰色+当前部分亮色） ====================
fprintf('--- 正在生成 %d 个子图进行分色可视化 ---\n', num_subplots);

% 定义颜色：每个子图的高亮部分用独特颜色，完整数据用灰色
highlight_colors = lines(num_subplots);  % 20种区分度高的亮色
full_data_color = [0.75 0.75 0.75];  % 完整数据的灰色（背景）

% 计算每份的点数（最后几份可能略多）
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

% 创建图形
figure('Color', 'w', 'Position', [100 100 1400 900]);  % 宽×高保持
sgtitle(sprintf('单个梯级分部分可视化 (%d - %d)', ...
                current_step_data.Start_loc(1), current_step_data.Start_loc(end)), ...
        'FontSize', 16, 'FontWeight', 'bold', 'Color', 'k');

% 循环绘制子图：每个子图 = 完整数据（灰色） + 当前部分（对应亮色）
for i = 1:num_subplots
    subplot(rows, cols, i);
    hold on; grid on;
    
    % 1. 绘制完整数据（灰色背景）
    scatter(current_step_data.Azimuth, current_step_data.Elevation, point_size, full_data_color, 'filled');
    
    % 2. 绘制当前子图对应的部分数据（亮色突出）
    current_part_data = current_step_data(part_indices{i}, :);
    scatter(current_part_data.Azimuth, current_part_data.Elevation, point_size, highlight_colors(i,:), 'filled');
    
    % 图形样式优化
    % 仅第一列子图显示y轴标签
    if mod(i-1, cols) == 0  % 第一列：i=1,6,11,16
        ylabel('仰角 (°)', 'Color', 'k', 'FontSize', 7);
    end
    % 仅最后一行子图显示x轴标签
    if i > rows*(cols-1)-1  % 最后一行：i=16~20
        xlabel('方位角 (°)', 'Color', 'k', 'FontSize', 7);
    end

    % 统一坐标轴范围和样式
    xlim([155 175])
    ylim([60 75])
    set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'Box', 'on', ...
        'FontSize', 6, 'TickDir', 'in');  % 保持样式一致
    axis equal;  % 保持坐标比例
    hold off;
end

disp('绘图完成!');

fprintf(' %.2f us , 每个图的分辨率为%.2f us \n', ...
    (current_step_data.Start_loc(end)-current_step_data.Start_loc(1))*5/1e3, ...
    (current_step_data.Start_loc(end)-current_step_data.Start_loc(1))*5/1e3/20);